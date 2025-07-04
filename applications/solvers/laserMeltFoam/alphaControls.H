const dictionary& alphaControls = mesh.solverDict(alphaM.name());

label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));
