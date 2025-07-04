/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

Class
    Foam::csvTableReader

Description
    Reads an interpolation table from a file - CSV-format

SourceFiles
    tableReader.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_csvTableReader_H
#define Foam_csvTableReader_H

#include "tableReader.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class csvTableReader Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class csvTableReader
:
    public tableReader<Type>
{
    // Private Data

        //- Does the file have a header line?
        const bool headerLine_;

        //- Column for reference (lookup) value
        const label refColumn_;

        //- Labels of the components
        const labelList componentColumns_;

        //- Separator character
        const char separator_;


    // Private Member Functions

        //- Get component columns entry
        static labelList getComponentColumns
        (
            const word& name,
            std::initializer_list<std::pair<const char*,int>> compat,
            const dictionary& dict
        );

        //- Read component values from the split string
        Type readValue(const UList<string>& strings) const;


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("csv");


    // Constructors

        //- Construct from dictionary
        csvTableReader(const dictionary& dict);

        //- Construct and return a copy
        virtual autoPtr<tableReader<Type>> clone() const
        {
            return autoPtr<tableReader<Type>>
            (
                new csvTableReader<Type>
                (
                    *this
                )
            );
        }


    //- Destructor
    virtual ~csvTableReader() = default;


    // Member Functions

        //- Read 1D table
        virtual void operator()
        (
            const fileName& fName,
            List<Tuple2<scalar, Type>>& data
        );

        //- Read 2D table - NotImplemented
        virtual void operator()
        (
            const fileName& fName,
            List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& data
        );

        //- Write the remaining parameters
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "csvTableReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
