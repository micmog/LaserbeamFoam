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

rho             1000;

    cp  800;
    cpsolid 600.0;
    kappa  30.0;
	kappasolid  20.0; 
	Tsolidus 1658;
	Tliquidus 1723;
    LatentHeat 2.7e5;
    beta    5.0e-6;


// ************************************************************************* //
