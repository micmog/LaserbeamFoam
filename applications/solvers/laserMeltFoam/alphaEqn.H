// {
//     // Temporarily making U relative to mesh motion
//     if (mesh.moving())
//     {
//         U -= fvc::reconstruct(mesh.phi());
//     }

//     // Updating alpha1
//     #include "alphaSuSp.H"
//     advector.advect(Sp, Su);

//     // Making U absolute again after advection step
//     if (mesh.moving())
//     {
//         U += fvc::reconstruct(mesh.phi());
//     }

//     #include "rhofs.H"
//     rhoPhi = advector.getRhoPhi(rho1f, rho2f);

//     alpha2 = 1.0 - alpha1;
//     mixture.correct();
// }

// scalar domainFraction = 0;
// if (porosityEnabled)
// {
//     const volScalarField& porosity = tporosity.cref();
//     rhoPhi *= scalar(1)/fvc::interpolate(porosity);
//     domainFraction = alpha1.weightedAverage(mesh.Vsc()*porosity).value();
// }
// else
// {
//     domainFraction = alpha1.weightedAverage(mesh.Vsc()).value();
// }

// Info<< "Phase-1 volume fraction = "
//     << domainFraction
//     << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
//     << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
//     << endl;

{
  // STEP 1: set up flux and velocity for Crank-Nicolson scheme

  // Initialize placeholders for current U and phi fields before the solution.
  // Note: needed for isoAdvector that uses truncated Crank-Nicolson scheme
  // (only flux and velocity fields are discretised with C-N scheme, while
  // alpha is always explicitly discretised).
  autoPtr<surfaceScalarField> curPhiPtr;
  autoPtr<volVectorField> curUPtr;

  // If topo changes have occured, we must not take the old flux
      Info<< "Setting up velocity and flux fields for isoAdvector." << endl;

      // Store current fields
      curPhiPtr.reset(new surfaceScalarField("tmpPhi", phi));
      curUPtr.reset(new volVectorField("tmpU", U));

      // Get Crank-Nicolson off-centre coefficient from isoAdvector subdict
      // Note: different definition as in ddtSchemes
      // - 0 = Euler explicit
      // - 0.5 = Full Crank-Nicolson
      // - 1 = Euler implicit
      const scalar CNCoeff = readScalar
      (
          mesh.solutionDict().subDict("isoAdvector").lookup
          (
              "CrankNicolsonCoeff"
          )
      );

      // Set current fields for Crank-Nicolson scheme
      phi = CNCoeff*phi + (1.0 - CNCoeff)*phi.oldTime();
      U = CNCoeff*U + (1.0 - CNCoeff)*U.oldTime();

  // STEP 2: Advect the free surface using isoAdvector

  // Create isoAdvection object
  isoAdvection advector(alphaM, phi, U);

  // Advect the free surface
  #include "alphaSuSp.H"
  advector.advect(Sp, Su);
      
  // STEP 3: Restore flux and velocity fields after isoAdvector
  if
  (
      curPhiPtr.valid()
   && curUPtr.valid()
  )
  {
      Info<< "Restoring velocity and flux fields after isoAdvector." << endl;

      phi = curPhiPtr();
      U = curUPtr();

      // Clear out stored fields
      curPhiPtr.clear();
      curUPtr.clear();
  }
  else if
  (
      (!curPhiPtr.valid() || !curUPtr.valid())
  )
  {
      // Sanity check
      FatalError<< "Current flux and velocity fields are invalid after the"
                << " isoAdvector advection step."
                << nl
                << "Something went wrong."
                << exit(FatalError);
  }

  Info<< "Material phase volume fraction = "
      << alphaM.weightedAverage(mesh.V()).value()
      << "  Min(alphaM) = " << min(alphaM).value()
      << "  Max(alphaM) = " << max(alphaM).value()
      << endl;

  interface.correct();

  // # include   "updateProperties.H"
  rho == alphaM*((gT*rhoL)+((scalar(1)-gT)*rhoS)) + (scalar(1) - alphaM)*rhoG;

  surfaceScalarField rhoMf (fvc::interpolate((gT*rhoL)+((scalar(1)-gT)*rhoS)));

  //Reverse steps to get phiAlpha
  surfaceScalarField phiAlpha ((advector.getRhoPhi(rhoMf, rhoG) - rhoG*phi)/(rhoMf-rhoG));
  
  rhoPhi = phiAlpha*(rhoMf-rhoG) + phi*rhoG;

  /*surfaceScalarField checkRhoPhi = phi*fvc::interpolate(alphaM)*(rhoMf-rhoG) + phi*rhoG;*/

  // rhoCpPhi = phiAlpha*(rhoMf*(fvc::interpolate(cM)) - rhoG*cpG) + phi*rhoG*cpG;
  // rhoCpPhi = fvc::interpolate((alphaM*cM)+((scalar(1)-alphaM)*cpG))*rhoPhi;
  // rhoMPhi = fvc::interpolate(alphaM)*rhoPhi * Hf;
  
  // Update interfaceFvDataVOF object after obtaining new alpha field
  // intFvDataVOF.update();
  
}
