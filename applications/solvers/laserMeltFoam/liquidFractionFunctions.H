/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Conection between volume fraction and temperature.

Author
    Filip Volaric, UCD

\*---------------------------------------------------------------------------*/
#ifndef liquidFractionFunctions_H
#define liquidFractionFunctions_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//temperature to liquid volume fraction - fields g = F(T)
tmp<volScalarField> liquidFraction
(
    const volScalarField& T,
    const dimensionedScalar& Tl,
    const dimensionedScalar& Ts
);

//volume fraction to the temperature - fields T = F^-1(g)
tmp<volScalarField> invLiquidFraction
(
    const volScalarField& gT,
    const dimensionedScalar& Tl,
    const dimensionedScalar& Ts
);

//derivative of liquid fraction by temperature dF/dT
tmp<volScalarField> DFbyDT
(
    const volScalarField& T,
    const dimensionedScalar& Tl,
    const dimensionedScalar& Ts
);

//[dF/dT]^-1
tmp<volScalarField> invDFbyDT
(
    const volScalarField& T,
    const dimensionedScalar& Tl,
    const dimensionedScalar& Ts
);

//temperature to liquid volume fraction - points g = F(T)
scalar liquidFraction(const scalar& T, const dimensionedScalar& Tl,
                const dimensionedScalar& Ts);

//liquid volume fraction to the temperature - points T = F^-1(g)
scalar invLiquidFraction(const scalar& gT, const dimensionedScalar& Tl,
                const dimensionedScalar& Ts);

//derivative of liquid fraction by temperature dF/dT
scalar DFbyDT(const scalar& T, const dimensionedScalar& Tl,
                const dimensionedScalar& Ts);
//[dF/dT]^-1
scalar invDFbyDT(const scalar& T, const dimensionedScalar& Tl,
                const dimensionedScalar& Ts);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} //end namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************** //
