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

nu              8e-07;

rho             2540;

	Tsolidus 843;
	Tliquidus 907;
    LatentHeat 3.98e5;
    beta    5.0e-6;
    poly_kappa   (155 0.01 0 0 0 0 0 0);
    poly_cp   (900.0 0.25 0 0 0 0 0 0);


// ************************************************************************* //
