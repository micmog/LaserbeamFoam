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

nu              5e-07;

rho             7500;

	Tsolidus 1658;
	Tliquidus 1723;
    LatentHeat 2.7e5;
    beta    5.0e-6;
    poly_kappa   (10 0.016 0 0 0 0 0 0);
    poly_cp   (450 0.16 0 0 0 0 0 0);


// ************************************************************************* //
