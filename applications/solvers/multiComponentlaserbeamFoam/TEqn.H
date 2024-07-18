{
    int iter = 0;

    dimensionedScalar residual("residual",dimensionSet(0, 0, 0, 0, 0), scalar(1.0));
    dimensionedScalar meanResidual("meanResidual",dimensionSet(0, 0, 0, 0, 0), scalar(1.0));
    dimensionedScalar residual2("residual2",dimensionSet(0, 0, 0, 0, 0), scalar(1.0));
    dimensionedScalar meanResidual2("meanResidual2",dimensionSet(0, 0, 0, 0, 0), scalar(1.0));

    TRHS *= 0.0;

    const volScalarField rhoCp(rho*cp);
    const surfaceScalarField rhophicp(fvc::interpolate(cp)*rhoPhi);


    do
    {
        iter++;
        epsilon1.storePrevIter();

        TRHS =
            LatentHeat
           *(
               fvc::ddt(rho,epsilon1)
             + fvc::div(rhoPhi, epsilon1)
           );

        TRHS.correctBoundaryConditions();

        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T)
          + fvm::div(rhophicp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
          - fvm::laplacian(kappa, T)
         ==
            fvModels.source(rhoCp,T)
          + Deposition
          - TRHS
          - Qv*mag(gradAlpha)
        );

        TEqn.relax();

        fvConstraints.constrain(TEqn);

        TEqn.solve();

        Tcorr = (TLiquidus-TSolidus)*epsilon1+TSolidus;

        epsilon1 = max(min(epsilon1+(epsilonRel*cp/LatentHeat)*(T-Tcorr), scalar(1)), scalar(0));

        residual = gMax(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField()));
        meanResidual = gSum(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField())*mesh.V())/gSum(mesh.V().field());

        Info<< "Correcting epsilon1, mean residual = " << meanResidual.value()
            << ", max residual = " << residual.value()
            << endl;
    }
    while ((iter < minTCorr || (residual > epsilonTol)) && iter <= maxTCorr);

    T.correctBoundaryConditions();
}
