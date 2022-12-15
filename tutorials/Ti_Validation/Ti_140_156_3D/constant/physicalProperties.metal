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

nu              2.35e-06;

rho             4420;


	Tsolidus 1880;
	Tliquidus 1928.15;
    LatentHeat 2.95653e5;
    beta    5.0e-4;
    
    poly_kappa   (5 0.012 0 0 0 0 0 0);
    poly_cp   (520 0.05 0 0 0 0 0 0);


// ************************************************************************* //
