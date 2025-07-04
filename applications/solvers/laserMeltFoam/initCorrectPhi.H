tmp<volScalarField> rAU;

if (correctPhi)
{
    rAU = new volScalarField
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rAU", dimTime/dimDensity, 1)
    );

    #include "correctPhi.H"
}
else
{
    CorrectPhi
    (
        U,
        phi,
        pd, //can we just switch p_rgh to pd here
        dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
        geometricZeroField(),
        pimple
    );

    #include "continuityErrs.H"
}
