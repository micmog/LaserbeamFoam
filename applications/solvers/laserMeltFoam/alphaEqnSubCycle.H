if (pimple.nCorrPIMPLE() > 1)
{
    // If nOuterCorrectors > 1 then for all but the first loop the advection
    // of alpha is done using an average, 0.5*phi+0.5*phiNew where phi is
    // the flux at the beginning of the time step and phiNew is the flux
    // estimate at the end of the time step from the previous outer
    // iteration. Similarly we use 0.5*U + 0.5*UNew in later iterations.
    if (pimple.firstIter())
    {
        // To recalculate the alpha1 update in subsequent iterations, we
        // must store its current value before overwriting with the new
        // value
        alphaM.storePrevIter();
        // Storing initial phi and U for use in later outer iterations.
        phi.storePrevIter();
        U.storePrevIter();
    }
    else
    {
        // Resetting alpha1 to value before advection in first PIMPLE
        // iteration.
        alphaM = alphaM.prevIter();

        // Temporarily setting U and phi with which to advect interface.

        // TODO all these U and phi fields need to be checked if they are correct, this was the old way:
        // U = 0.5*U.prevIter() + 0.5*U;
        // phi = 0.5*phi.prevIter() + 0.5*phi;


        // TODO should it be like this now? This is now in alphaEqn.H, no need here?
        // const scalar CNCoeff = readScalar
        // (
        //     mesh.solutionDict().subDict("isoAdvector").lookup
        //     (
        //         "CrankNicolsonCoeff"
        //     )
        // );
        // Set current fields for Crank-Nicolson scheme
        // phi = CNCoeff*phi + (1.0 - CNCoeff)*phi.oldTime();
        // U = CNCoeff*U + (1.0 - CNCoeff)*U.oldTime();

    }
}

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        mesh.newIOobject("rhoPhiSum"),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alphaM, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

if (!pimple.firstIter())
{

    // Resetting U and phi to value at latest iteration. 

    // TODO should this be different now?
    U = 2.0*U - U.prevIter();
    phi = 2.0*phi - phi.prevIter();

    // const scalar CNCoeff = readScalar
    // (
    //     mesh.solutionDict().subDict("isoAdvector").lookup
    //     (
    //         "CrankNicolsonCoeff"
    //     )
    // );

    // phi = CNCoeff*phi + (1.0 - CNCoeff)*phi.oldTime();
    // U = CNCoeff*U + (1.0 - CNCoeff)*U.oldTime();
}


// Info<< "Material phase volume fraction = "
// << alphaM.weightedAverage(mesh.V()).value()
// << "  Min(alphaM) = " << min(alphaM).value()
// << "  Max(alphaM) = " << max(alphaM).value()
// << endl;

// interface.correct();


// rho == alphaM*((gT*rhoL)+((scalar(1)-gT)*rhoS)) + (scalar(1) - alphaM)*rhoG;
// surfaceScalarField rhoMf = fvc::interpolate((gT*rhoL)+((scalar(1)-gT)*rhoS));

// //Reverse steps to get phiAlpha
// surfaceScalarField phiAlpha = (advector.getRhoPhi(rhoMf, rhoG) - rhoG*phi)/(rhoMf-rhoG);

// rhoPhi = phiAlpha*(rhoMf-rhoG) + phi*rhoG;

// // rho == alpha1*rho1 + alpha2*rho2;
