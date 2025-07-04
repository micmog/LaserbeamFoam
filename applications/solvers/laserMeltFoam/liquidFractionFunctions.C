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

Author
    Filip Volaric, UCD Dublin

Modified by:
    Gowthaman Parivendhan, UCD

\*---------------------------------------------------------------------------*/

#include "liquidFractionFunctions.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  tmp<volScalarField> liquidFraction(
      const volScalarField &T,
      const dimensionedScalar &Tl,
      const dimensionedScalar &Ts)
  {
    tmp<volScalarField> tresult(
        new volScalarField(
            IOobject(
                "F(" + T.name() + ")",
                T.time().timeName(),
                T.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            T.mesh(),
            dimensionedScalar("0", dimless, 0)));

    volScalarField &result = tresult.ref();
    
    scalarField &resultI = result.internalFieldRef();

    const scalarField &TI = T.internalField();

    forAll(resultI, cellI)
    {
      resultI[cellI] = liquidFraction(TI[cellI], Tl, Ts);
    }

    forAll(result.boundaryField(), patchI)
    {
      if (
          T.boundaryField()[patchI].type() != emptyFvPatchField<vector>::typeName)
      {
        scalarField &pResult = result.boundaryFieldRef()[patchI];
        const scalarField &pT = T.boundaryField()[patchI];

        forAll(pResult, faceI)
        {
          pResult[faceI] = liquidFraction(pT[faceI], Tl, Ts);
        }
      }
    }

    // result.correctBoundaryConditions();

    return tresult;
  }

  tmp<volScalarField> invLiquidFraction(
      const volScalarField &gT,
      const dimensionedScalar &Tl,
      const dimensionedScalar &Ts)
  {
    tmp<volScalarField> tresult(
        new volScalarField(
            IOobject(
                "invF(" + gT.name() + ")",
                gT.time().timeName(),
                gT.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            gT.mesh(),
            dimensionedScalar("0", dimTemperature, 0)));

    volScalarField &result = tresult.ref();
    scalarField &resultI = result.internalFieldRef();

    const scalarField &gTI = gT.internalField();

    forAll(resultI, cellI)
    {
      resultI[cellI] = invLiquidFraction(gTI[cellI], Tl, Ts);
    }

    forAll(result.boundaryField(), patchI)
    {
      if (
          gT.boundaryField()[patchI].type() != emptyFvPatchField<vector>::typeName)
      {
        scalarField &pResult = result.boundaryFieldRef()[patchI];
        const scalarField &pgT = gT.boundaryField()[patchI];

        forAll(pResult, faceI)
        {
          pResult[faceI] = invLiquidFraction(pgT[faceI], Tl, Ts);
        }
      }
    }

    // result.correctBoundaryConditions();

    return tresult;
  }

  tmp<volScalarField> DFbyDT(
      const volScalarField &T,
      const dimensionedScalar &Tl,
      const dimensionedScalar &Ts)
  {
    tmp<volScalarField> tresult(
        new volScalarField(
            IOobject(
                "dF(" + T.name() + ")/dT",
                T.time().timeName(),
                T.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            T.mesh(),
            dimensionedScalar("0", dimless / dimTemperature, 0)));

    volScalarField &result = tresult.ref();
    scalarField &resultI = result.internalFieldRef();

    const scalarField &TI = T.internalField();

    forAll(resultI, cellI)
    {
      resultI[cellI] = DFbyDT(TI[cellI], Tl, Ts);
    }

    forAll(result.boundaryField(), patchI)
    {
      if (
          T.boundaryField()[patchI].type() != emptyFvPatchField<vector>::typeName)
      {
        scalarField &pResult = result.boundaryFieldRef()[patchI];
        const scalarField &pT = T.boundaryField()[patchI];

        forAll(pResult, faceI)
        {
          pResult[faceI] = DFbyDT(pT[faceI], Tl, Ts);
        }
      }
    }

    // result.correctBoundaryConditions();

    return tresult;
  }

  tmp<volScalarField> invDFbyDT(
      const volScalarField &T,
      const dimensionedScalar &Tl,
      const dimensionedScalar &Ts)
  {
    tmp<volScalarField> tresult(
        new volScalarField(
            IOobject(
                "dF(" + T.name() + ")/dT^-1",
                T.time().timeName(),
                T.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            T.mesh(),
            dimensionedScalar("0", dimTemperature / dimless, 0)));

    volScalarField &result = tresult.ref();
    scalarField &resultI = result.internalFieldRef();

    const scalarField &TI = T.internalField();

    forAll(resultI, cellI)
    {
      resultI[cellI] = invDFbyDT(TI[cellI], Tl, Ts);
    }

    forAll(result.boundaryField(), patchI)
    {
      if (
          T.boundaryField()[patchI].type() != emptyFvPatchField<vector>::typeName)
      {
        scalarField &pResult = result.boundaryFieldRef()[patchI];
        const scalarField &pT = T.boundaryField()[patchI];

        forAll(pResult, faceI)
        {
          pResult[faceI] = invDFbyDT(pT[faceI], Tl, Ts);
        }
      }
    }

    // result.correctBoundaryConditions();

    return tresult;
  }

  scalar liquidFraction(
      const scalar &T,
      const dimensionedScalar &Tl,
      const dimensionedScalar &Ts)
  {
    scalar result = 0;

    if (T >= Ts.value() && T <= Tl.value())
    {
      result = (T - Ts.value()) / (Tl.value() - Ts.value() + SMALL);
    }
    else if (T > Tl.value())
    {
      result = 1.0;
    }
    return result;
  }

  scalar invLiquidFraction(
      const scalar &gT,
      const dimensionedScalar &Tl,
      const dimensionedScalar &Ts)
  {
    scalar result = Ts.value();

    if (gT > 0.0 && gT < 1.0)
    {
      result = Ts.value() + gT * (Tl.value() - Ts.value());
    }
    else if (gT >= 1.0)
    {
      result = Tl.value();
    }

    return result;
  }

  scalar DFbyDT(
      const scalar &T,
      const dimensionedScalar &Tl,
      const dimensionedScalar &Ts)
  {
    scalar result = 0;

    if (T >= Ts.value() && T <= Tl.value())
    {
      result = scalar(1.0) / (Tl.value() - Ts.value() + SMALL);
    }

    if (result > Foam::pow(10.0, 10))
    {
      result = Foam::pow(10.0, 10);
    }

    return result;
  }

  scalar invDFbyDT(
      const scalar &T,
      const dimensionedScalar &Tl,
      const dimensionedScalar &Ts)
  {
    scalar result = 0 + SMALL;

    if (T >= Ts.value() && T <= Tl.value())
    {
      result = Tl.value() - Ts.value();
    }

    return result;
  }
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
