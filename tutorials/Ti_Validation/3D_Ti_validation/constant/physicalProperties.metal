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

    cp  750;
    cpsolid 600.0;
    kappa  16.0;
	kappasolid  8.0; 
	Tsolidus 1890;
	Tliquidus 1928.15;
    LatentHeat 2.95653e5;
    beta    5.0e-6;


// ************************************************************************* //
