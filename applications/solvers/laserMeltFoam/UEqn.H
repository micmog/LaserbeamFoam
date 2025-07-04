surfaceScalarField muEff(
    "muEff",
    threePhaseProperties.muf() + fvc::interpolate(rho * turbulence->nut()));

pVap = min(max((T-(TVap-(TSmooth/2.0)))/TSmooth,0.0),1.0)*0.54 * p0 * Foam::exp(Hv*Mm*((T-TVap)/(R*T*TVap)));

dimensionedScalar deltaN = 1e-8/pow(average(mesh.V()), 1.0/3.0);
volVectorField N (fvc::grad(alphaM)/(mag(fvc::grad(alphaM)) + deltaN));
volScalarField damper (2.0 * rho /(rhoM+rhoG) ); //To account for high density difference
A = -1.0 * C * alphaM * (Foam::pow(scalar(1)-gT,2))/(Foam::pow(gT,3)+SMALL);

MRF.correctBoundaryVelocity(U);

// fvVectorMatrix UEqn( // old eq
//     fvm::ddt(rho, U) + fvm::div(rhoPhi, U) + MRF.DDt(rho, U) + turbulence->divDevRhoReff(rho, U) ==
//     fvOptions(rho, U));

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(muEff, U)
  - (fvc::grad(U) & fvc::grad(muEff))  - fvm::Sp(A,U)
  - (DSigbyDT * (fvc::grad(T) - N*(N & fvc::grad(T)))) * mag(fvc::grad(alphaM)) * damper == (rhoT-rho)*g
  //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
  //== -Foam::max(alphaM * gT * betaM * (T-Tl),0.0) * g * rhoM
);

// #include "UEqnAddPorosity.H"

UEqn.relax();

fvOptions.constrain(UEqn);

// if (pimple.momentumPredictor())
// {
//   solve(
//       UEqn ==
//       fvc::reconstruct(
//           (
//               mixture.surfaceTensionForce() - ghf * fvc::snGrad(rho) - fvc::snGrad(p_rgh)) *
//           mesh.magSf()));

//   fvOptions.correct(U);
// }


if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
             fvc::interpolate((interface.sigmaK()+pVap)*damper)*fvc::snGrad(alphaM)
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(pd)
             )*mesh.magSf()
        )
    );
}
