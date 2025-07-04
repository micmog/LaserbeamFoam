{
  // scalar relRes = GREAT; // not used anywhere
  scalar perfRes = GREAT;

  SolverPerformance<scalar> Tp;

  em = (gT * emL) + ((1.0 - gT) * emS);
  rhoCp = rho * ((alphaM * cM) + ((1.0 - alphaM) * cpG));

  rhoCp.correctBoundaryConditions();

  rhoCpPhi = fvc::interpolate((alphaM * cM) + ((1.0 - alphaM) * cpG)) * rhoPhi;

  rhoMPhi = fvc::interpolate(alphaM) * rhoPhi * Hf;

  volScalarField delGradAlpha(mag(fvc::grad(alphaM)) * 2.0 * (rhoCp) / (rhoM * cM + rhoG * cpG));
  // Qv = 0.82*Hv*Mm*p0*Foam::exp(Hv*Mm*((T-TVap)/(R*T*TVap)))/Foam::pow(2*M_PI*Mm*R*T,0.5);
  // volScalarField fluxVc = (-hC*(T-TRef) - SB*em*(Foam::pow(T,4) - Foam::pow(TRef,4)) - Qv)*delGradAlpha;

  for (int i = 0; i < nTCorrectors; i++)
  {

    Info << "Temperature solver Iteration :" << i + 1 << endl;
    T.storePrevIter();
    gT.storePrevIter();

    #include "source.H"
    // volScalarField DQvbyDT = (Hv*Mm/(R*Foam::pow(T,2))) - (1/(2*T));

    Qv = min(max((T - (TVap - (TSmooth / 2.0))) / TSmooth, 0.0), 1.0) * 0.82 * Hv * Mm * p0 * Foam::exp(Hv * Mm * ((T - TVap) / (R * T * TVap))) / Foam::pow(2 * M_PI * Mm * R * T, 0.5);

    volScalarField fluxVc ((hC * TRef + SB * em * (Foam::pow(TRef, 4) + 3 * Foam::pow(T, 4)) - Qv) * delGradAlpha);
    // volScalarField fluxVc = (-hC*(T-TRef) - SB*em*(Foam::pow(T,4) - Foam::pow(TRef,4)) - Qv)*delGradAlpha;
    volScalarField fluxCc (-(hC + 4 * SB * em * Foam::pow(T, 3)) * delGradAlpha);
    
    if (ray_tracing_on) // TODO  - this could be better
    { 
       laser_source = laser.deposition();
    }
    else
    {
       laser_source = laserFlux/delta;
    }

    fvScalarMatrix TEqn(fvm::ddt(rhoCp, T)  + fvm::div(rhoCpPhi, T, "div(phi,T)")  - fvm::SuSp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)  - fvm::Sp(SI, T) - fvm::Sp(fluxCc, T) + fvc::div(rhoMPhi, gT) - fvm::laplacian(kEff, T) == SE + fluxVc + laser_source);

    TEqn.relax();

    Tp = TEqn.solve();
    #include "liquidFraction.H"

    // perfRes = Tp.initialResidual();
    scalar absResidualT(gMax(mag(T.internalField() - T.prevIter().internalField()).ref()));
    scalar gTResidual(gMax(mag(gT.internalField() - gT.prevIter().internalField()).ref()));

    Info << "T Absolute residual: " << absResidualT << endl;
    Info << "gT Absolute residual: " << gTResidual << endl;

    if ((gTResidual < gTTol && perfRes < temperatureTol) || (gTResidual < gTTol && absResidualT < absTol))
    {
      break;
    }
    else if (i == (nTCorrectors - 1))
    {
      Warning 
          << "Temperature equation did not converge!" << endl;
    }

  }

  // Update gT after TEqn update
  gT = liquidFraction(T, Tl, Ts);
  alphaL = alphaM * gT;
  ddtGT = fvc::ddt(gT);
  T.correctBoundaryConditions();

  // #include "updateProperties.H"
  // T.correctBoundaryConditions();

  Info << "min/max(T) = "
       << gMin(T) << ", " << gMax(T) << endl;
}
