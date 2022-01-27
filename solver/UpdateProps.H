    cp=(epsilon1*(alpha1*cp1 + alpha2*cp2)) + ((1.0-epsilon1)*(alpha1*cp1solid + alpha2*cp2solid));
    kappa=(epsilon1*(alpha1*kappa1 + alpha2*kappa2)) + ((1.0-epsilon1)*(alpha1*kappa1solid + alpha2*kappa2solid));
    TSolidus=alpha1*Tsolidus1 + alpha2*Tsolidus2;
    TLiquidus=alpha1*Tliquidus1 + alpha2*Tliquidus2;
    LatentHeat=alpha1*LatentHeat1 + alpha2*LatentHeat2;
    beta=alpha1*beta1 + alpha2*beta2;

    rhok=(1.0-max(epsilon1*(beta)*(Temperature-TSolidus),0.0));
    rhok.correctBoundaryConditions();



    Dummyalpha=fvc::average(alpha1);

    volVectorField gradepsilon1=fvc::grad(Dummyalpha);
    gradepsilon1.correctBoundaryConditions();

    nneps1=fvc::grad(Dummyalpha)/((mag(gradepsilon1) + deltaN));
    nneps1.correctBoundaryConditions();


    volVectorField gradalpha1=fvc::grad(alpha1);
    gradalpha1.correctBoundaryConditions();

    n_alpha=fvc::grad(alpha1)/((mag(gradalpha1) + deltaN));
    n_alpha.correctBoundaryConditions();


/*
forAll( mesh.C(), celli)
{
if(alpha1[celli]>1e-5){
TSolidus[celli]=Tsolidus1.value();
TLiquidus[celli]=Tliquidus1.value();
}
else{
TSolidus[celli]=alpha1[celli]*Tsolidus1.value() + alpha2[celli]*Tsolidus2.value();
TLiquidus[celli]=alpha1[celli]*Tliquidus1.value() + alpha2[celli]*Tliquidus2.value();
}
}


forAll( mesh.C(), celli)
{
if(alpha1[celli]<0.99){
LatentHeat[celli]=LatentHeat2.value();

}
else{
LatentHeat[celli]=LatentHeat1.value();
}
}
*/




volScalarField e1temp = fvc::average(epsilon1);



forAll( mesh.C(), celli)
{
if(e1temp[celli]<=0.95){
epsilon1mask[celli]=0.0;
}
else{
epsilon1mask[celli]=epsilon1[celli];
}

}

epsilon1mask.correctBoundaryConditions();


    gh=epsilon1mask*(g & mesh.C());
    ghf= fvc::interpolate(epsilon1mask)*(g & mesh.Cf());
