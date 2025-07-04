/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "interpolationTable.H"
#include "openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::interpolationTable<Type>::readTable()
{
    // preserve the original (unexpanded) fileName to avoid absolute paths
    // appearing subsequently in the write() method
    fileName fName(fileName_);

    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorInFunction
            << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are okay
    check();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTable<Type>::interpolationTable()
:
    List<value_type>(),
    bounding_(bounds::repeatableBounding::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(nullptr)
{}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable
(
    const List<Tuple2<scalar, Type>>& values,
    const bounds::repeatableBounding bounding,
    const fileName& fName
)
:
    List<value_type>(values),
    bounding_(bounding),
    fileName_(fName),
    reader_(nullptr)
{}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable(const fileName& fName)
:
    List<value_type>(),
    bounding_(bounds::repeatableBounding::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>())
{
    readTable();
}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable(const dictionary& dict)
:
    List<value_type>(),
    bounding_
    (
        bounds::repeatableBoundingNames.getOrDefault
        (
            "outOfBounds",
            dict,
            bounds::repeatableBounding::WARN,
            true  // Failsafe behaviour
        )
    ),
    fileName_(dict.get<fileName>("file")),
    reader_(tableReader<Type>::New(dict))
{
    readTable();
}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable
(
     const interpolationTable& tbl
)
:
    List<value_type>(tbl),
    bounding_(tbl.bounding_),
    fileName_(tbl.fileName_),
    reader_(tbl.reader_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::interpolationTable<Type>::check() const
{
    const List<value_type>& list = *this;

    scalar prevValue(0);

    label i = 0;
    for (const auto& item : list)
    {
        const scalar& currValue = item.first();

        // Avoid duplicate values (divide-by-zero error)
        if (i && currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
        ++i;
    }
}


template<class Type>
void Foam::interpolationTable<Type>::write(Ostream& os) const
{
    os.writeEntry("file", fileName_);
    os.writeEntry("outOfBounds", bounds::repeatableBoundingNames[bounding_]);
    if (reader_)
    {
        reader_->write(os);
    }
}


template<class Type>
Type Foam::interpolationTable<Type>::rateOfChange(scalar lookupValue) const
{
    const List<value_type>& list = *this;

    const label n = list.size();

    if (n <= 1)
    {
        // Not enough entries for a rate of change
        return Zero;
    }

    const scalar minLimit = list.first().first();
    const scalar maxLimit = list.last().first();

    if (lookupValue < minLimit)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")\n"
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")\n"
                    << "    Zero rate of change." << endl;

                // Behaviour as per CLAMP
                return Zero;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                return Zero;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // Adjust lookupValue to >= minLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue - minLimit, span) + minLimit;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")\n"
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")\n"
                    << "    Zero rate of change." << endl;

                // Behaviour as per CLAMP
                return Zero;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                return Zero;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // Adjust lookupValue <= maxLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue - minLimit, span) + minLimit;
                break;
            }
        }
    }

    label lo = 0;
    label hi = 0;

    // Look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= list[i].first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        return Zero;
    }
    else if (hi == 0)
    {
        // This treatment should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            (list[hi].second() - list[lo].second())
          / (list[hi].first() + minLimit - list[lo].first())
        );
    }


    // Normal rate of change
    return
    (
        (list[hi].second() - list[lo].second())
      / (list[hi].first() - list[lo].first())
    );
}


template<class Type>
Type Foam::interpolationTable<Type>::interpolateValue
(
    const List<Tuple2<scalar, Type>>& list,
    scalar lookupValue,
    bounds::repeatableBounding bounding
)
{
    const label n = list.size();

    if (n <= 1)
    {
        #ifdef FULLDEBUG
        if (!n)
        {
            FatalErrorInFunction
                << "Cannot interpolate from zero-sized table" << nl
                << exit(FatalError);
        }
        #endif

        return list.first().second();
    }

    const scalar minLimit = list.first().first();
    const scalar maxLimit = list.last().first();

    if (lookupValue < minLimit)
    {
        switch (bounding)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")\n"
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")\n"
                    << "    Continuing with the first entry" << endl;

                // Behaviour as per CLAMP
                return list.first().second();
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                return list.first().second();
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // adjust lookupValue to >= minLimit
                const scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue - minLimit, span) + minLimit;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (bounding)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")\n"
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")\n"
                    << "    Continuing with the last entry" << endl;

                // Behaviour as per 'CLAMP'
                return list.last().second();
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                return list.last().second();
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // Adjust lookupValue <= maxLimit
                const scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue - minLimit, span) + minLimit;
                break;
            }
        }
    }


    label lo = 0;
    label hi = 0;

    // Look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= list[i].first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        return list[hi].second();
    }
    else if (hi == 0)
    {
        // This treatment should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            list[lo].second()
          + (list[hi].second() - list[lo].second())
          * (lookupValue / minLimit)
        );
    }


    // Normal interpolation
    return
    (
        list[lo].second()
      + (list[hi].second() - list[lo].second())
      * (lookupValue - list[lo].first())
      / (list[hi].first() - list[lo].first())
    );
}


template<class Type>
Type Foam::interpolationTable<Type>::interpolateValue
(
    scalar lookupValue
) const
{
    return interpolateValue(*this, lookupValue, bounding_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::interpolationTable<Type>::interpolateValues
(
    const UList<scalar>& vals
) const
{
    auto tfld = tmp<Field<Type>>::New(vals.size());
    auto& fld = tfld.ref();

    forAll(fld, i)
    {
        fld[i] = interpolateValue(vals[i]);
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::interpolationTable<Type>::operator=
(
    const interpolationTable<Type>& rhs
)
{
    if (this == &rhs)
    {
        return;
    }

    static_cast<List<value_type>&>(*this) = rhs;
    bounding_ = rhs.bounding_;
    fileName_ = rhs.fileName_;
    reader_.reset(rhs.reader_.clone());
}


template<class Type>
const Foam::Tuple2<Foam::scalar, Type>&
Foam::interpolationTable<Type>::operator[](label idx) const
{
    const List<value_type>& list = *this;
    const label n = list.size();

    if (n <= 1)
    {
        idx = 0;

        #ifdef FULLDEBUG
        if (!n)
        {
            FatalErrorInFunction
                << "Cannot interpolate from zero-sized table" << nl
                << exit(FatalError);
        }
        #endif
    }
    else if (idx < 0)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "index (" << idx << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "index (" << idx << ") underflow" << nl
                    << "    Continuing with the first entry" << nl;

                // Behaviour as per 'CLAMP'
                idx = 0;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                idx = 0;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                while (idx < 0)
                {
                    idx += n;
                }
                break;
            }
        }
    }
    else if (idx >= n)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "index (" << idx << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "index (" << idx << ") overflow" << nl
                    << "    Continuing with the last entry" << nl;

                // Behaviour as per 'CLAMP'
                idx = n - 1;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                idx = n - 1;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                while (idx >= n)
                {
                    idx -= n;
                }
                break;
            }
        }
    }

    return list[idx];
}


template<class Type>
Type Foam::interpolationTable<Type>::operator()(scalar lookupValue) const
{
    return interpolateValue(*this, lookupValue, bounding_);
}


// ************************************************************************* //
