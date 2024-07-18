
// cp=(epsilon1*(alpha1*cp1 + (1.0-alpha1)*cp2)) + ((1.0-epsilon1)*(alpha1*cp1solid + (1.0-alpha1)*cp2solid));
// kappa=(epsilon1*(alpha1*kappa1 + (1.0-alpha1)*kappa2)) + ((1.0-epsilon1)*(alpha1*kappa1solid + (1.0-alpha1)*kappa2solid));
TSolidus=alpha1*Tsolidus1 + (1.0-alpha1)*Tsolidus2;
TLiquidus=alpha1*Tliquidus1 + (1.0-alpha1)*Tliquidus2;
// LatentHeat=alpha1*LatentHeat1 + (1.0-alpha1)*LatentHeat2;
beta=alpha1*beta1 + (1.0-alpha1)*beta2;

rhok=(1.0-max(epsilon1*(beta)*(T-TSolidus),0.0));
rhok.correctBoundaryConditions();

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





//update the cell size fields
const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

const vectorField& CI = mesh.C();
scalarField& yDimI = laser.yDim();

forAll(CI, celli)
{
    // vector XYZ = mesh.C()[celli];
    // xcoord[celli]=XYZ.x();
    // zcoord[celli]=XYZ.z();

    const cell& cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll(pLabels, pointi)
    {
        pLocal[pointi] = pp[pLabels[pointi]];
    }

    // xDim[celli] = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
    yDimI[celli] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
    // zDim[celli] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
}

// xDim.correctBoundaryConditions();
laser.yDim().correctBoundaryConditions();
// zDim.correctBoundaryConditions();
//update the cell size fields










// refineflag=fvc::average(alpha1);
metalaverage = fvc::average(alpha1);

laser.refineFlag() *= 0.0;
laser.refineFlag() += metalaverage;

// forAll( mesh.C(), celli)
// {
// if((metalaverage[celli]>0.1&&metalaverage[celli]<0.9)||refineflag2[celli]>0.1){
// refineflag2[celli]=0.5;
// }
// }

alpha_filtered=alpha1;

forAll(CI, celli)
{
    if (alpha1[celli] > 0.99)
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

forAll(CI, celli)
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
forAll(CI, celli)
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

gh = epsilon1mask*(g & mesh.C());
ghf = fvc::interpolate(epsilon1mask)*(g & mesh.Cf());
