    IOdictionary solutionDict
    (
        IOobject
        (
        "fvSolution",
        runTime.system(), //folder where the fvSolution file is located 
        runTime,
        IOobject::MUST_READ
        )
    );



    dictionary& meltingDict = solutionDict.subDict("MELTING");

    int minTCorr(readLabel(meltingDict.lookup("minTempCorrector")));
    int maxTCorr(readLabel(meltingDict.lookup("maxTempCorrector")));

    dimensionedScalar epsilonTol(readScalar(meltingDict.lookup("epsilonTolerance")));
    dimensionedScalar epsilonRel(readScalar(meltingDict.lookup("epsilonRelaxation")));

    bool damperSwitch = meltingDict.lookupOrDefault<bool>("damperSwitch", false);
