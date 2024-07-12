TSolidus = mixture.Tsolidus();
TLiquidus = mixture.Tliquidus();
Tvapourisation = mixture.Tvap();

kappa = (epsilon1*mixture.kappa())+((1.0-epsilon1)*mixture.kappasolid());
cp = (epsilon1*mixture.Cp())+((1.0-epsilon1)*mixture.Cpsolid());

molarmass = mixture.molarmass();
resistivity = mixture.elec_resistivity();
e_density = mixture.e_num_density();

beta = mixture.beta();
LatentHeat = mixture.Latentheat();

LatentHeatGas = mixture.Latentheatgas();

rhok = (1.0 - max(epsilon1*(mixture.beta())*(T - (mixture.Tsolidus())), 0.0));

condensate *= 0.0;

forAllConstIter(PtrDictionary<phase>, mixture.phases(), iter)
{
    const phase& phase1 = iter();
    const volScalarField& alpha1 = phase1;

    // Info<<"name, "<<" , "<<alpha1.name()<<endl;
    if (phase1.name() != "air")
    {
        condensate += alpha1;
    }
}

condensate = max(min(condensate,1.0),0.0);

alpha_filtered = condensate;

forAll(mesh.C(), celli)
{
    // dont care about optical properties of the gas phase for absorptivity
    resistivity[celli] = max(resistivity[celli],0.0)/(condensate[celli]+deltaN.value()); 
    e_density[celli] = max(e_density[celli],0.0)/(condensate[celli]+deltaN.value());
    Tvapourisation[celli] = Tvapourisation[celli]/(condensate[celli]+deltaN.value());

    if (condensate[celli] > 0.99)
    {
        alpha_filtered[celli] = 1.0;
    }
    else
    {
        if (condensate[celli] < 0.01)
        {
            alpha_filtered[celli] = 0.0;
        }
        else
        {
            alpha_filtered[celli] = condensate[celli];
        }

    }
}
alpha_filtered.correctBoundaryConditions();





const volVectorField gradAlpha(fvc::grad(condensate));

// pVap = 0.54 * p0 * Foam::exp(LatentHeatGas*molarmass*((T-Tvapourisation)/(constant::physicoChemical::R*T*Tvapourisation)));

// Qv = 0.82*LatentHeatGas*molarmass*p0*Foam::exp(LatentHeatGas*molarmass
//     *((T-Tvapourisation)/(constant::physicoChemical::R*T*Tvapourisation)))/Foam::pow(2*constant::mathematical::pi*molarmass*constant::physicoChemical::R*T,0.5);


pVap *= 0.0;
Qv *= 0.0;
forAll( mesh.C(), celli)
{

if(T[celli]>Tvapourisation[celli]){
pVap[celli] =  0.54 * p0.value() * Foam::exp(LatentHeatGas[celli]*molarmass[celli]*((T[celli]-Tvapourisation[celli])/(constant::physicoChemical::R.value()*T[celli]*Tvapourisation[celli])));

Qv[celli] = 0.82*LatentHeatGas[celli]*molarmass[celli]*p0.value()*Foam::exp(LatentHeatGas[celli]*molarmass[celli]
    *((T[celli]-Tvapourisation[celli])/(constant::physicoChemical::R.value()*T[celli]*Tvapourisation[celli])))/Foam::pow(2*constant::mathematical::pi*molarmass[celli]*constant::physicoChemical::R.value()*T[celli],0.5);

}

// if(condensate[celli]<0.5){//Dont want recoil effects in gas phase, just at interface
// pVap[celli] *= 0.0;
// Qv[celli] *= 0.0;
// }
}



volVectorField gradepsilon1(fvc::grad(alpha_filtered));

n_filtered = fvc::grad(alpha_filtered)/((mag(gradepsilon1) + deltaN));


const volScalarField e1temp(fvc::average(epsilon1));


forAll(mesh.C(), celli)
{
    if (e1temp[celli] <= 0.95)
    {
        epsilon1mask[celli] = 0.0;
    }
    else
    {
        epsilon1mask[celli] = epsilon1[celli];
    }
}

epsilon1mask.correctBoundaryConditions();

gh = epsilon1*(g & mesh.C());
ghf =  fvc::interpolate(epsilon1)*(g & mesh.Cf());
