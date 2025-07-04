#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

//- Overset specific

// Add solver-specific interpolations
{
    wordHashSet& nonInt =
        const_cast<wordHashSet&>(Stencil::New(mesh).nonInterpolatedFields());

    nonInt.insert("HbyA");
    nonInt.insert("grad(p_rgh)");
    nonInt.insert("nHat");
    nonInt.insert("surfaceIntegrate(nHatf)");
    nonInt.insert("surfaceIntegrate(phi+meshPhi)");
    nonInt.insert("surfaceIntegrate(phi)");
    nonInt.insert("surfaceIntegrate(phiHbyA)");
    nonInt.insert("surfaceSum(((S|magSf)*S)");
    nonInt.insert("surfaceIntegrate(((rAUf*magSf)*snGradCorr(p_rgh)))");
    nonInt.insert("cellMask");
    nonInt.insert("cellDisplacement");
    nonInt.insert("interpolatedCells");
    nonInt.insert("cellInterpolationWeight");
    nonInt.insert("pcorr");
}


// Mask field for zeroing out contributions on hole cells
#include "createCellMask.H"

// Create bool field with interpolated cells
#include "createInterpolatedCells.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
