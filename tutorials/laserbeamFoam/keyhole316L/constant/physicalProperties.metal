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
    object      physicalProperties.metal;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

viscosityModel  constant;

//nu              6.67e-04;
nu              4.8e-07; // 4.8e-07

rho             7950;

Tsolidus		1658;
Tliquidus		1727;
LatentHeat		245e3;

beta			4.95e-5;
    
poly_kappa  	(18  0 0 0 0 0 0 0);
poly_cp			(600 0 0 0 0 0 0 0);


// ************************************************************************* //
