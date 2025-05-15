{
    int iter = 0;
    scalarField residual(epsilon1.primitiveField().size(), 1.0);
    scalar maxResidual = 1.0;
    scalar meanResidual = 1.0;
    scalar Tresidual = 1.0;

    TRHS *= 0.0;

    const volScalarField rhoCp(rho*cp);
    const surfaceScalarField rhophicp(fvc::interpolate(cp)*rhoPhi);

    // const volTensorField gradU(fvc::grad(U));
    // volTensorField tau
    // (
    //     ((mixture.mu()*gradU) + (mixture.mu()*gradU.T()))
    //   - ((2.0/3.0)*mixture.mu()*(fvc::div(phi))*I)
    // );

    do
    {
        epsilon1.storePrevIter();
        T.storePrevIter();
        TRHS =
            LatentHeat*((fvc::ddt(rho, epsilon1) + fvc::div(rhoPhi, epsilon1)));

        Qv =
            0.82*LatentHeatVap*Mm*p0
           *Foam::exp(LatentHeatVap*Mm*((T - Tvap)/(R*T*Tvap)))
           /Foam::pow(2*M_PI*Mm*R*T,0.5);

        if (damperSwitch)
        {
            forAll(mesh.C(), celli)
            {
                thermalDamper[celli] =
                    2.0*rhoCp[celli]
                   /(
                       (rho1.value()*polycp_m.value(T[celli]))
                     + (rho2.value()*polycp_g.value(T[celli]))
                   );
            }

            thermalDamper.correctBoundaryConditions();
        }

        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T)
          + fvm::div(rhophicp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
          - fvm::laplacian(kappa, T)
          //- (tau && gradU) // Viscous dissipation
         ==
            fvModels.source(rhoCp,T)
          + laser.deposition()
          - Qv*mag(gradAlpha)*thermalDamper
          - TRHS
        );

        TEqn.relax();

        fvConstraints.constrain(TEqn);

        TEqn.solve();

        T.relax();

        Tcorr = (TLiquidus - TSolidus)*epsilon1 + TSolidus;

        epsilon1 =
            max
            (
                min
                (
                    epsilon1 + (epsilonRel*cp/LatentHeat)*(T - Tcorr),
                    scalar(1)
                ),
                scalar(0)
            );

        residual =
            mag
            (
                epsilon1.primitiveField()
              - epsilon1.prevIter().primitiveField()
            );

        maxResidual = gMax(residual);

        meanResidual = gSum(residual*mesh.V())/gSum(mesh.V().field());

        Info<< "Correcting epsilon1, mean residual = " << meanResidual
            << ", max residual = " << maxResidual
            << endl;

        ddte1 = fvc::ddt(epsilon1);
    }
    while
    (
        (
            ++iter < minTCorr
         // || (meanResidual > epsilonTol.value())
         || (maxResidual > epsilonTol.value())
        )
     && iter <= maxTCorr
    );

    gradT = fvc::grad(T);
}





