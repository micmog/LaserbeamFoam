
{


    int iter = 0;

    dimensionedScalar residual("residual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar meanResidual("meanResidual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar Tresidual("Tresidual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));

    TRHS*=0.0;

    volScalarField rhoCp = rho*cp;
    rhoCp.correctBoundaryConditions();
    surfaceScalarField rhoPhiCp=fvc::interpolate(cp)*rhoPhi;
    surfaceScalarField rhophicp(fvc::interpolate(cp)*rhoPhi);

        do{
            iter++;
            epsilon1.storePrevIter();
            T.storePrevIter();
            TRHS = LatentHeat*((fvc::ddt(rho,epsilon1) + fvc::div(rhoPhi, epsilon1)));
            TRHS.correctBoundaryConditions();

            volTensorField gradU = fvc::grad(U);
            volTensorField tau = ((mixture.mu()*gradU) + (mixture.mu()*gradU.T()))-((2.0/3.0)*mixture.mu()*(fvc::div(phi))*I);

            Qv = min(max((T-(Tvap-(TSmooth/2.0)))/TSmooth,0.0),1.0)*0.82*LatentHeatVap*Mm*p0*Foam::exp(LatentHeatVap*Mm*((T-Tvap)/(R*T*Tvap)))/Foam::pow(2*M_PI*Mm*R*T,0.5);
            Qv.correctBoundaryConditions();

            if(damperSwitch)
            {
                thermalDamper = 2.0*(rhoCp)/(rho1*cp1 + rho2*cp2);
            }
   
            // ViscousDissipation=(tau && gradU);

            fvScalarMatrix  TEqn
            (
            fvm::ddt(rhoCp,T)
            + fvm::div(rhophicp, T)
            - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
            - fvm::laplacian(kappa, T)
       //     -(tau && gradU) //Viscous dissipation
            ==
            fvModels.source(rhoCp,T) + (Deposition) - Qv*mag(gradAlpha)*thermalDamper - (TRHS)
            );

            TEqn.relax();
            fvConstraints.constrain(TEqn);

            TEqn.solve();

        // fvConstraints.correct(T);

        Tcorr=(TLiquidus-TSolidus)*epsilon1+TSolidus;
        Tcorr.correctBoundaryConditions();

        epsilon1=max(min(epsilon1+(epsilonRel*cp/LatentHeat)*(T-Tcorr),scalar(1)),scalar(0));
        epsilon1.correctBoundaryConditions();

        residual = gMax(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField()));

        meanResidual = gSum(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField())*mesh.V())/gSum(mesh.V().field());


        Info << "Correcting epsilon1, mean residual = " << meanResidual.value()
         << ", max residual = " << residual.value()
	     << endl;

        ddte1=fvc::ddt(epsilon1);
        ddte1.correctBoundaryConditions();
        }while((iter < minTCorr || (residual > epsilonTol) ) && iter <= maxTCorr );

    T.correctBoundaryConditions();
    gradT=fvc::grad(T);
    gradT.correctBoundaryConditions();
    }





