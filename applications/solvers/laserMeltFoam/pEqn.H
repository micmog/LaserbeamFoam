// {
//     if (correctPhi)
//     {
//         rAU.ref() = 1.0/UEqn.A();
//     }
//     else
//     {
//         rAU = 1.0/UEqn.A();
//     }

//     surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
//     volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
//     surfaceScalarField phiHbyA
//     (
//         "phiHbyA",
//         fvc::flux(HbyA)
//       + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
//     );
//     MRF.makeRelative(phiHbyA);

//     if (p_rgh.needReference())
//     {
//         fvc::makeRelative(phiHbyA, U);
//         adjustPhi(phiHbyA, U, p_rgh);
//         fvc::makeAbsolute(phiHbyA, U);
//     }

//     surfaceScalarField phig
//     (
//         (
//             mixture.surfaceTensionForce()
//           - ghf*fvc::snGrad(rho)
//         )*rAUf*mesh.magSf()
//     );

//     phiHbyA += phig;

//     // Update the pressure BCs to ensure flux consistency
//     constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

//     while (pimple.correctNonOrthogonal())
//     {
//         fvScalarMatrix p_rghEqn
//         (
//             fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
//         );

//         p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
//         p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

//         if (pimple.finalNonOrthogonalIter())
//         {
//             phi = phiHbyA - p_rghEqn.flux();

//             p_rgh.relax();

//             U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
//             U.correctBoundaryConditions();
//             fvOptions.correct(U);
//         }
//     }

//     #include "continuityErrs.H"

// Correct Uf if the mesh is moving
// fvc::correctUf(Uf, U, phi);

//     // Make the fluxes relative to the mesh motion
//     fvc::makeRelative(phi, U);

//     p == p_rgh + rho*gh;

//     if (p_rgh.needReference())
//     {
//         p += dimensionedScalar
//         (
//             "p",
//             p.dimensions(),
//             pRefValue - getRefCellValue(p, pRefCell)
//         );
//         p_rgh = p - rho*gh;
//     }

//     if (!correctPhi)
//     {
//         rAU.clear();
//     }
// }

{
  volScalarField rUA(1.0 / UEqn.A());
  surfaceScalarField rUAf(fvc::interpolate(rUA));

  U = rUA * UEqn.H();

  surfaceScalarField phiU(
      "phiU",
      (fvc::interpolate(U) & mesh.Sf())
      // + fvc::ddtPhiCorr(rUA, rho, U, phi)
      // + fvc::fvcDdtPhiCorr(U, phi) // TODO this needs to be added somehow
  );

  adjustPhi(phiU, U, pd);

  phi = phiU +
        (fvc::interpolate((interface.sigmaK() + pVap) * damper) * fvc::snGrad(alphaM)
         //+ fvc::interpolate(recoil)
         //+ fvc::interpolate((pVap*N)*mag(fvc::grad(alphaM))*damper)
         - ghf * fvc::snGrad(rho)) *
            rUAf * mesh.magSf();

  while (pimple.correctNonOrthogonal())
  {
    pd.storePrevIter();

    // volScalarField rhoSource = alphaM*(fvc::ddt(rhoM)+(fvc::grad(rhoM)&U))/rhoM /*+ (1-alphaM)*fvc::ddt(rhoG)/rhoG*/;
    fvScalarMatrix pdEqn(
        fvm::laplacian(rUAf, pd) == fvc::div(phi) /*+ rhoSource*/
    );

    // pdEqn.setReference(pdRefCell, pdRefValue);
    pdEqn.setReference(pdRefCell, pdRefValue);

    pdEqn.relax();

    pdEqn.solve(
        // mesh.solutionDict().solution(pd.select(pimple.finalInnerIter()))
        mesh.solutionDict(pd.select(pimple.finalInnerIter())) // TODO this might not work right
    );

    if (pimple.finalNonOrthogonalIter())
    {
      phi -= pdEqn.flux();
      pd.relax();
      // U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
      // U.correctBoundaryConditions();
    }
  }

  U += rUA * fvc::reconstruct((phi - phiU) / rUAf);
  U.correctBoundaryConditions();
}
