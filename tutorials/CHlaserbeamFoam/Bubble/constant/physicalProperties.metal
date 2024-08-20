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

nu              1e-06;

//mu	9.97e-4;

rho             997.0;


    sigma_E   (3.46e6 0.0 0 0 0 0 0 0);
    mu_M   (0.8e-6 0.0 0 0 0 0 0 0);
    
    
	Tsolidus 1658;
	Tliquidus 1723;
    LatentHeat 2.7e5;
    beta    5.0e-6;
    poly_kappa   (25 0.0 0 0 0 0 0 0);
    poly_cp   (700 0.0 0 0 0 0 0 0);


// ************************************************************************* //
