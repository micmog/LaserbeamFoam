IOdictionary laserProperties(
    IOobject(
        "laserProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

const Switch ray_tracing_on = laserProperties.lookupOrDefault<Switch>("ray_tracing_on", false);

if(ray_tracing_on==true){
  Info << "\nRay-tracing: ON" <<endl;
  Info << "Using Ray-tracing laser model" <<endl;
  Info << "Don't forget to orient the mesh correctly (laser must be in the positive y axis)!\n" <<endl;
}else{
  Info << "\nRay-tracing: OFF" <<endl;
  Info << "Using default laser model\n" <<endl;
}

interpolationTable<scalar> laserPowerSeries(laserProperties.subDict("timeVsLaserPower"));
interpolationTable<vector> laserPositionSeries(laserProperties.subDict("timeVsLaserPosition"));

const scalar kappaM(readScalar(laserProperties.lookup("kappaM"))); // Absorptivity
const scalar emL(readScalar(laserProperties.lookup("emL")));
const scalar emS(readScalar(laserProperties.lookup("emS")));

volScalarField em(
    IOobject(
        "em",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0));

// const scalar eps(readScalar(laserProperties.lookup("eps"))); //Emissivity

word laserPatchName(word(laserProperties.lookup("laserPatchName")));

const dimensionedScalar r0("laserRadius",  dimensionSet(0, 1, 0, 0, 0, 0, 0) ,laserProperties.lookup("laserRadius")); // Laser radius
dimensionedScalar t("t", dimensionSet(0, 0, 1, 0, 0, 0, 0), 0.0);
dimensionedScalar Power("Power", dimensionSet(1, 2, -3, 0, 0, 0, 0), 0.0);
dimensionedVector laserPos("laserPos", dimensionSet(0, 1, 0, 0, 0, 0, 0), vector(0.0, 0.0, 0.0));

volScalarField px(mesh.C().component(0)); // Pointer for this, unnecessary, stays the same throughout simulation
volScalarField py(mesh.C().component(1));

volScalarField flb(
    IOobject(
        "flb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0));

volScalarField delta(
    IOobject(
        "delta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0),
    "zeroGradient");

dimensionedScalar SB(
    "SB",
    dimensionSet(1, 0, -3, -4, 0, 0, 0),
    scalar(5.67e-8)); // Stefan-Boltzmann Constant

volScalarField laserFlux(
    IOobject(
        "laserFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("0", dimensionSet(1, 0, -3, 0, 0, 0, 0), 0.0));

    