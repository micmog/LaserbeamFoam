    MRF.correctBoundaryVelocity(U);


dimensionedScalar DarcyConstantlarge("DarcyConstantlarge",dimensionSet(1,-3,-1,0,0),scalar(1.0e6));//1e10
dimensionedScalar DarcyConstantsmall("DarcyConstantsmall",dimensionSet(0,0,0,0,0),scalar(1.0e-12));//1e-6




if(PowderSim){

DC=(DarcyConstantlarge*Foam::pow((1.0-epsilon1mask),2)/(Foam::pow(epsilon1mask,3)+DarcyConstantsmall));
DC.correctBoundaryConditions();
}
else{
    DC=(DarcyConstantlarge*Foam::pow((1.0-epsilon1),2)/(Foam::pow(epsilon1,3)+DarcyConstantsmall));
    DC.correctBoundaryConditions();
}

Marangoni = epsilon1mask*mixture.MarangoniForce(T);


    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevTau(rho, U)
      + fvm::Sp(DC, U)
      + Marangoni
     ==
        fvModels.source(rho, U)
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                    +fvc::interpolate(pVap)*fvc::snGrad(condensate)
                  - ghf*fvc::snGrad(rho*rhok)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvConstraints.constrain(U);
    }
