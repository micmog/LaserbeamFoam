/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "laserThreePhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

  //- Calculate and return the laminar viscosity

  // TODO this works somehow but needs to be fixed
  word threePhaseMixture::getPhaseName(const word &key) const
  {

    // ERROR: no member named 'isDict' in 'Foam::word'

    // if (key.isDict())
    // {
    //   return key;
    // }
    // else
    // {
    //   return word.lookup(key);
    // }
    return key;
  }

  void threePhaseMixture::calcNu()
  {
    nuModel1_->correct();
    nuModel2_->correct();
    volScalarField limitedAlphaM(
        "limitedAlphaM",
        min(max(alphaM_, scalar(0)), scalar(1)));
    volScalarField limitedgT(
        "limitedgT",
        min(max(gT_, scalar(0)), scalar(1)));
    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu() / (limitedAlphaM * ((limitedgT * rhoL_) + ((scalar(1) - limitedgT) * rhoS_)) + (scalar(1) - limitedAlphaM) * rhoG_);
  }

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  threePhaseMixture::threePhaseMixture(
        const volVectorField &U,
        const surfaceScalarField &phi,
        const word &alphaMName,
        const word &gTName,
        const word &TName):
        transportModelNew(U, phi),
        phase1Name_(getPhaseName("material")),
        phase2Name_(getPhaseName("gas")),
        nuModel1_(
            viscosityModel::New(
                "nu1",
                subDict(phase1Name_),
                U,
                phi)),
        nuModel2_(
            viscosityModel::New(
                "nu2",
                subDict(phase2Name_),
                U,
                phi)),
        rhoL_("rhoL",nuModel1_->viscosityProperties()),
        rhoS_("rhoS",nuModel1_->viscosityProperties()),
        rhoG_("rhoG",nuModel2_->viscosityProperties()),
        U_(U),
        phi_(phi),
        alphaM_(U_.db().lookupObject<const volScalarField>(alphaMName)),
        gT_(U_.db().lookupObject<const volScalarField>(gTName)),
        T_(U_.db().lookupObject<const volScalarField>(TName)),
        Ts_("Ts", dimensionSet(0, 0, 0, 1, 0, 0, 0), 0.0),
        Tl_("Tl", dimensionSet(0, 0, 0, 1, 0, 0, 0), 0.0),
        nu_(
            IOobject(
                "nu",
                U_.time().timeName(),
                U_.db()),
            U_.mesh(),
            dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
            calculatedFvPatchScalarField::typeName)
  {
    IOdictionary thermalProperties(
        IOobject(
            "thermalProperties",
            U_.time().constant(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE));
            dimensionedScalar Ts_ ("Ts", thermalProperties.subDict("metalProperties"));
            dimensionedScalar Tl_ ("Tl", thermalProperties.subDict("metalProperties"));
            calcNu();
  }

  // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
  tmp<volScalarField> threePhaseMixture::alphaM() const
  {
    return tmp<volScalarField>(
        new volScalarField(
            "alpha.material",
            alphaM_));
  }

  tmp<volScalarField> threePhaseMixture::mu() const
  {
    volScalarField limitedAlphaM(min(max(alphaM_, scalar(0)), scalar(1)));
    volScalarField limitedgT(min(max(gT_, scalar(0)), scalar(1)));

    return tmp<volScalarField>(
        new volScalarField(
            "mu_threePhaseMixture",
            limitedAlphaM * ((limitedgT * rhoL_) + ((scalar(1) - limitedgT) * rhoS_)) * nuModel1_->nu() + (scalar(1) - limitedAlphaM) * rhoG_ * nuModel2_->nu()));
  }

  tmp<surfaceScalarField> threePhaseMixture::muf() const
  {
    surfaceScalarField alphaMf(min(max(fvc::interpolate(alphaM_), scalar(0)), scalar(1)));

    surfaceScalarField gTf(min(max(fvc::interpolate(gT_), scalar(0)), scalar(1)));

    return tmp<surfaceScalarField>(
        new surfaceScalarField(
            "muf_threePhaseMixture",
            alphaMf * ((gTf * rhoL_) + ((scalar(1) - gTf) * rhoS_)) * fvc::interpolate(nuModel1_->nu()) + (scalar(1) - alphaMf) * rhoG_ * fvc::interpolate(nuModel2_->nu())));
  }

  tmp<surfaceScalarField> threePhaseMixture::nuf() const
  {
    surfaceScalarField alphaMf(min(max(fvc::interpolate(alphaM_), scalar(0)), scalar(1)));

    surfaceScalarField gTf(min(max(fvc::interpolate(gT_), scalar(0)), scalar(1)));

    return tmp<surfaceScalarField>(
        new surfaceScalarField(
            "nuf_threePhaseMixture",
            (
                alphaMf * ((gTf * rhoL_) + ((scalar(1) - gTf) * rhoS_)) * fvc::interpolate(nuModel1_->nu()) + (scalar(1) - alphaMf) * rhoG_ * fvc::interpolate(nuModel2_->nu())) /
                (alphaMf * ((gTf * rhoL_) + ((scalar(1) - gTf) * rhoS_)) + (scalar(1) - alphaMf) * rhoG_)));
  }

  bool threePhaseMixture::read()
  {
    if (transportModelNew::read())
    {
      if (
          nuModel1_().read(subDict(phase1Name_)) &&
          nuModel2_().read(subDict(phase2Name_)))
      {
        nuModel1_->viscosityProperties().readEntry("rhoS", rhoS_);
        nuModel1_->viscosityProperties().readEntry("rhoL", rhoL_);
        nuModel2_->viscosityProperties().readEntry("rhoG", rhoG_);
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
