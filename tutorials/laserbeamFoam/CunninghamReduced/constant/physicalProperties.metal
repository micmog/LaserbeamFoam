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

nu              3.9e-06;

rho             4457; // Ti-alloy

Tsolidus		1878; // Ti-alloy
Tliquidus		1928; // Ti-alloy
LatentHeat		2.9e5; // Ti-alloy

beta			4.95e-5; // Ti-alloy
    
poly_kappa 	 	(30 0 0 0 0 0 0 0); // Ti-alloy
poly_cp			(700 0 0 0 0 0 0 0); // Ti-alloy


// ************************************************************************* //
