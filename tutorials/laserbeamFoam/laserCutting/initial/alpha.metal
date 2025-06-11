/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  3.0.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    object      alpha.metal.orig;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

	#includeEtc "caseDicts/setConstraintTypes"

    leftWall
    {
		type            zeroGradient;
		/*
		type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
        */
    }
    
    rightWall
    {
		type            zeroGradient;
		/*
		type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
        */
    }

	frontAndBack
	{
        type            zeroGradient;
        /*
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
        */
	}
	
    lowerWall
    {
		type            zeroGradient;
    }
    
    atmosphere
    {
    	type            zeroGradient;
    	/*
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
        */
    }
     
}


// ************************************************************************* //
