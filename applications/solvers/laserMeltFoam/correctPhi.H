CorrectPhi
(
    U,
    phi,
    pd, // can we just switch p_rgh to pd here
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"
