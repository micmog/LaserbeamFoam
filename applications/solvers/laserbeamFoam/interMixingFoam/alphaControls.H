const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaCorr(alphaControls.get<label>("nAlphaCorr"));

label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));
