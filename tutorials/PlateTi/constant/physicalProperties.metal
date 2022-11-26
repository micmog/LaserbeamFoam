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

nu              4.0e-06;

rho             4420;

    cp  1100;
    cpsolid 650.0;
    kappa  40.0;
	kappasolid  30.0; 
	Tsolidus 1890;
	Tliquidus 1941.15;
    LatentHeat 2.95653e5;
    beta    7.0e-5;


// ************************************************************************* //
