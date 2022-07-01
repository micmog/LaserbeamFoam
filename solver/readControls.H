    #include "readTimeControls.H"

    const dictionary& meltingDict = mesh.solutionDict().subDict("MELTING");
    int minTCorr(readLabel(meltingDict.lookup("minTempCorrector")));
    int maxTCorr(readLabel(meltingDict.lookup("maxTempCorrector")));

    dimensionedScalar epsilonTol(readScalar(meltingDict.lookup("epsilonTolerance")));
    dimensionedScalar epsilonRel(readScalar(meltingDict.lookup("epsilonRelaxation")));

    bool damperSwitch = meltingDict.lookupOrDefault<bool>("damperSwitch", false);



IOdictionary LaserProperties
(
    IOobject
    (
        "LaserProperties",    // dictionary name
        runTime.constant(),     // dict is found in "constant"
        mesh,                   // registry for the dict
        IOobject::MUST_READ,    // must exist, otherwise failure
        IOobject::NO_WRITE      // dict is only read by the solver
    )
);



    scalar HS_a(readScalar(LaserProperties.lookup("HS_a")));
    scalar HS_bg(readScalar(LaserProperties.lookup("HS_bg")));
    scalar HS_velocity(readScalar(LaserProperties.lookup("HS_velocity")));
    scalar HS_lg(readScalar(LaserProperties.lookup("HS_lg")));
    scalar HS_Q(readScalar(LaserProperties.lookup("HS_Q")));
    vector V_incident(LaserProperties.lookup("V_incident"));

    scalar wavelength(readScalar(LaserProperties.lookup("wavelength")));
    scalar e_num_density(readScalar(LaserProperties.lookup("e_num_density")));
    scalar elec_resistivity(readScalar(LaserProperties.lookup("elec_resistivity"))); // this is temperature dependent - will include this in future versions
      

    

    // scalar HS_deposition_cutoff(readScalar(meltingDict.lookup("HS_deposition_cutoff")));
    // scalar Oscillation_Amplitude(readScalar(meltingDict.lookup("Oscillation_Amplitude")));
    // scalar Oscillation_Frequency(readScalar(meltingDict.lookup("Oscillation_Frequency")));


