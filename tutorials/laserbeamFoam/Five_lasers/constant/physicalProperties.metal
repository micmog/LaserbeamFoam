/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

viscosityModel  constant;

nu              1.35e-06;

rho             4420;


	Tsolidus 1878;
	Tliquidus 1928.15;
    LatentHeat 2.95653e5;
    beta    5.0e-4;
    
    poly_kappa   (10 0.015 0 0 0 0 0 0);
    //poly_cp   (244.8 9.587e-1 -3.77e-4 6.5e-8 -4.14e-12 0 0 0);
    poly_cp   (520 0.075 0 0 0 0 0 0);
    elec_resistivity	2.0e-6;

// ************************************************************************* //
