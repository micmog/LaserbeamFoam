
cp=(epsilon1*(alpha1*cp1 + (1.0-alpha1)*cp2)) + ((1.0-epsilon1)*(alpha1*cp1solid + (1.0-alpha1)*cp2solid));
kappa=(epsilon1*(alpha1*kappa1 + (1.0-alpha1)*kappa2)) + ((1.0-epsilon1)*(alpha1*kappa1solid + (1.0-alpha1)*kappa2solid));
TSolidus=alpha1*Tsolidus1 + (1.0-alpha1)*Tsolidus2;
TLiquidus=alpha1*Tliquidus1 + (1.0-alpha1)*Tliquidus2;
LatentHeat=alpha1*LatentHeat1 + (1.0-alpha1)*LatentHeat2;
beta=alpha1*beta1 + (1.0-alpha1)*beta2;

rhok=(1.0-max(epsilon1*(beta)*(T-TSolidus),0.0));
rhok.correctBoundaryConditions();



alpha_filtered=alpha1;

forAll( mesh.C(), celli)
{
    if(alpha1[celli]>0.99)
    {

        alpha_filtered[celli]=1.0;

    }
    else
    {

        if(alpha1[celli]<0.01)
        {

            alpha_filtered[celli]=0.0;

        }
        else
        {
            alpha_filtered[celli]=alpha1[celli];
        }

    }
}
alpha_filtered.correctBoundaryConditions();

const volVectorField gradepsilon1(fvc::grad(alpha_filtered));

n_filtered = fvc::grad(alpha_filtered)/((mag(gradepsilon1) + deltaN));

const volScalarField e1temp(fvc::average(epsilon1));



forAll( mesh.C(), celli)
{
    if(e1temp[celli]<=0.95)
    {
        epsilon1mask[celli]=0.0;
    }
    else
    {
        epsilon1mask[celli]=epsilon1[celli];
    }

}

epsilon1mask.correctBoundaryConditions();


gh=epsilon1mask*(g & mesh.C());
ghf= fvc::interpolate(epsilon1mask)*(g & mesh.Cf());
