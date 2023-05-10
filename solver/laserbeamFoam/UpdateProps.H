
// cp=(epsilon1*(alpha1*cp1 + (1.0-alpha1)*cp2)) + ((1.0-epsilon1)*(alpha1*cp1solid + (1.0-alpha1)*cp2solid));
// kappa=(epsilon1*(alpha1*kappa1 + (1.0-alpha1)*kappa2)) + ((1.0-epsilon1)*(alpha1*kappa1solid + (1.0-alpha1)*kappa2solid));
TSolidus=alpha1*Tsolidus1 + (1.0-alpha1)*Tsolidus2;
TLiquidus=alpha1*Tliquidus1 + (1.0-alpha1)*Tliquidus2;
// LatentHeat=alpha1*LatentHeat1 + (1.0-alpha1)*LatentHeat2;
beta=alpha1*beta1 + (1.0-alpha1)*beta2;

rhok=(1.0-max(epsilon1*(beta)*(T-TSolidus),0.0));
rhok.correctBoundaryConditions();


const volScalarField rhoCp(rho*cp);


        if (damperSwitch)
        {
        damper = 2.0 * rho /(rho1+rho2); //To account for high density difference

        forAll( mesh.C(), celli)
{

thermalDamper[celli]=2.0*rhoCp[celli]/((rho1.value()*polycp_m.value(T[celli]))+(rho2.value()*polycp_g.value(T[celli])));


}


        }


forAll( mesh.C(), celli)
{
if(alpha1[celli]>0.05){
TSolidus[celli]=Tsolidus1.value();
TLiquidus[celli]=Tliquidus1.value();
}
}

forAll( mesh.C(), celli)
{
kappa[celli]=(max(min(alpha1[celli],1.0),0.0)*(polykappa_m.value(T[celli])))+((1.0-max(min(alpha1[celli],1.0),0.0))*polykappa_g.value(T[celli]));
cp[celli]=(max(min(alpha1[celli],1.0),0.0)*(polycp_m.value(T[celli])))+((1.0-max(min(alpha1[celli],1.0),0.0))*polycp_g.value(T[celli]));
}

kappa.correctBoundaryConditions();


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

LatentHeat=alpha_filtered*LatentHeat1 + (1.0-alpha_filtered)*LatentHeat2;

volVectorField gradepsilon1(fvc::grad(alpha_filtered));

gradepsilon1.correctBoundaryConditions();

n_filtered = fvc::grad(alpha_filtered)/((mag(gradepsilon1) + deltaN));

n_filtered.correctBoundaryConditions();

// volScalarField 
epsilon1.correctBoundaryConditions();
volScalarField e1temp(fvc::average(epsilon1));

e1temp.correctBoundaryConditions();

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

volScalarField e2temp(fvc::average(epsilon1mask));
epsilon1mask2*=0.0;
forAll( mesh.C(), celli)
{
    if(e2temp[celli]<=0.95)
    {
        epsilon1mask2[celli]=0.0;
    }
    else
    {
        epsilon1mask2[celli]=epsilon1mask[celli];
    }

}



const volVectorField gradAlpha = (fvc::grad(alpha1));

const volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN));
Marangoni = epsilon1*Marangoni_Constant*(gradT-nHatM*(nHatM & gradT))*mag(gradAlpha);




gh=epsilon1mask*(g & mesh.C());
ghf= fvc::interpolate(epsilon1mask)*(g & mesh.Cf());
