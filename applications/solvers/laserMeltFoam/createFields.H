#include "createRDeltaT.H"

Info << "Reading field pd\n"
     << endl;
volScalarField pd(
    IOobject(
        "pd",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

#include "readThermalProperties.H"

volScalarField cpS(
    IOobject(
        "cpS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 2, -2, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

Info << "Reading field alpha\n"
     << endl;

volScalarField alphaM(
    IOobject(
        "alpha.material",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

Info << "Reading field T\n"
     << endl;
volScalarField T(
    IOobject(
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

volScalarField gT(
    IOobject(
        "gT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    liquidFraction(T, Tl, Ts),
    zeroGradientFvPatchScalarField::typeName);

Info << "Reading field alpha.liquid\n"
     << endl;
volScalarField alphaL(
    IOobject(
        "alpha.liquid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    alphaM *gT,
    zeroGradientFvPatchScalarField::typeName);

Info << "After alphaL update:\t" << "  Min(alphaL) = " << min(alphaL).value()
     << "  Max(alphaL) = " << max(alphaL).value() << endl;

Info << "After gT update:\t" << "  Min(gT) = " << min(gT).value()
     << "  Max(gT) = " << max(gT).value() << endl;

Info << "Reading field U\n"
     << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

#include "createPhi.H"

Info << "Reading transportProperties\n"
     << endl;

threePhaseMixture threePhaseProperties(U, phi, "alpha.material", "gT", "T");

volScalarField rhoM(
    IOobject(
        "rhoM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0, 0, 0), 7950.0),
    zeroGradientFvPatchScalarField::typeName);

volScalarField cM(
    IOobject(
        "cM",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 2, -2, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

volScalarField kM(
    IOobject(
        "kM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ),
    mesh,
    dimensionedScalar("0", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

volScalarField kEff(
    IOobject(
        "kEff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT),
    alphaM *(kM) + (scalar(1) - alphaM) * kG,
    zeroGradientFvPatchScalarField::typeName);

const dimensionedScalar &rhoL = threePhaseProperties.rhoL();
const dimensionedScalar &rhoS = threePhaseProperties.rhoS();
const dimensionedScalar &rhoG = threePhaseProperties.rhoG();

// Need to store rho for ddt(rho, U)
volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE),
    alphaM *((gT * rhoL) + ((scalar(1) - gT) * rhoS)) + (scalar(1) - alphaM) * rhoG,
    zeroGradientFvPatchScalarField::typeName);

rho.oldTime();

volScalarField rhoCp(
    IOobject(
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT),
    rho *(alphaM *cM + (scalar(1) - alphaM) * cpG),
    zeroGradientFvPatchScalarField::typeName);

volScalarField A(
    IOobject(
        "A",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT),
    -1.0 * C * alphaM * (Foam::pow(scalar(1) - gT, 2)) / (Foam::pow(gT, 3) + SMALL),
    zeroGradientFvPatchScalarField::typeName);

volScalarField rhoT(
    IOobject(
        "rhoT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

#include "updateProperties.H"

Info << "max/min(rho): " << max(rho).value() << ", " << min(rho).value() << endl;

Info << "max/min(cpS): " << max(cpS).value() << ", " << min(cpS).value() << endl;
Info << "max/min(cM): " << max(cM).value() << ", " << min(cM).value() << endl;

Info << "max/min(kM): " << max(kM).value() << ", " << min(kM).value() << endl;
Info << "max/min(kEff): " << max(kEff).value() << ", " << min(kEff).value() << endl;

volScalarField deltaH(
    IOobject(
        "deltaH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    rhoM *(cpL - cpS) * T + rhoM * Hf);

Info << "Creating field implicit source Sp\n"
     << endl;
volScalarField SI(
    IOobject(
        "SI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(1, -1, -3, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

Info << "Creating field explicit source Sc\n"
     << endl;
volScalarField SE(
    IOobject(
        "SE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

  volScalarField laser_source(
    IOobject(
        "laser_source",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

dictionary controls = mesh.solutionDict().subDict("Controls");

// ------------------------- ------------------------- ------------------------- ------------------------- -------------------------
// TODO All these variables are used in either TEq.H or liquidFraction.H but Clang gives a warning that they are not used anywhere
//      so this [[maybe_unused]] is just to suppress the warnings

[[maybe_unused]] double relax = controls.lookupOrDefault<double>("relaxCorr", 1.0);
[[maybe_unused]] label nTCorrectors(readLabel(pimple.dict().lookup("nTCorrectors")));
[[maybe_unused]] scalar temperatureTol(readScalar(pimple.dict().subDict("residualControl").subDict("temperatureTol").lookup("tolerance")));
[[maybe_unused]] scalar absTol(readScalar(pimple.dict().subDict("residualControl").subDict("absTol").lookup("tolerance")));
[[maybe_unused]] scalar gTTol(readScalar(pimple.dict().subDict("residualControl").subDict("gTTol").lookup("tolerance")));
// ------------------------- ------------------------- ------------------------- ------------------------- -------------------------

surfaceScalarField rhoCpPhi(
    IOobject(
        "rhoCpPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    fvc::interpolate(rhoCp) * phi);

surfaceScalarField rhoMPhi(
    IOobject(
        "rhoMPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    fvc::interpolate(alphaM *((gT * rhoL) + ((scalar(1) - gT) * rhoS))) * phi * Hf);

Info << "Calculating field g.h\n"
     << endl;

label pdRefCell = 0;
scalar pdRefValue = 0.0;

setRefCell(pd, pimple.dict(), pdRefCell, pdRefValue);

// Mass flux
surfaceScalarField rhoPhi(
    IOobject(
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    fvc::interpolate(rho) * phi);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p(
    IOobject(
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    pd + rho * gh);

scalar pRefValue = 0.0;

if (pd.needReference())
{
  pRefValue = readScalar(pimple.dict().lookup("pRefValue"));

  p += dimensionedScalar(
      "p",
      p.dimensions(),
      pRefValue - getRefCellValue(p, pdRefCell));
}

mesh.setFluxRequired(pd.name());
mesh.setFluxRequired(alphaM.name());

#include "createMRF.H"
#include "createFvOptions.H"
#include "createPorosity.H"

interfaceProperties interface(alphaM, U, threePhaseProperties);

const dimensionedScalar DSigbyDT("STgrad", threePhaseProperties);

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModelNew> turbulence(incompressible::turbulenceModelNew::New(U, phi, threePhaseProperties));

volScalarField pVap(
    IOobject(
        "pVap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName); // Recoil Pressure

volScalarField Qv(
    IOobject(
        "Qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(1, 0, -3, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

// Residual Stress calculation
volVectorField gradTSol(
    IOobject(
        "gradTSol",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("zero", dimless, vector::zero));

volScalarField solidificationTime(
    IOobject(
        "solidificationTime",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("minusOne", dimless, -1),
    zeroGradientFvPatchScalarField::typeName);

volScalarField ddtGT(fvc::ddt(gT));

// Adding ray tracing laser

laserHeatSource laser(mesh);

volScalarField alpha_filtered
(
    IOobject
    (
        "alpha_filtered",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_filtered",dimensionSet(0, 0, 0, -0, 0),0.0)
);

volVectorField n_filtered
(
    IOobject
    (
        "n_filtered",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("n_filtered",dimensionSet(0, 0, 0, 0, 0),vector::zero)
);


const dimensionedScalar elec_resistivity_metal ("elec_resistivity", metal); //Molar mass of metal

// dimensionedScalar elec_resistivity_metal("elec_resistivity",dimensionSet(1,3,-3,0,0),transportProperties_metal.lookup("elec_resistivity"));

        volScalarField electrical_resistivity//electrical resistivity of condensate only
    (
        IOobject
        (
            "electrical_resistivity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("electrical_resistivity",dimensionSet(1,3,-3,0,0),0.0)
    );

electrical_resistivity = elec_resistivity_metal;

volScalarField metalaverage
(
    IOobject
    (
        "metalaverage",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("metalaverage",dimensionSet(0,0,0,0,0),0.0)
);

dimensionedScalar deltaN=1e-8/pow(average(mesh.V()), 1.0/3.0);

metalaverage = fvc::average(alphaM);

