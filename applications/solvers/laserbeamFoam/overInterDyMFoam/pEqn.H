{
    rAU = 1.0/UEqn.A();

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField H("H", UEqn.H());

    volVectorField HbyA("HbyA", U);

    HbyA = constrainHbyA(rAU*H, U, p_rgh);

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    MRF.makeRelative(phiHbyA);

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(cellMask*rho)
        )*rAUf*faceMask*mesh.magSf()
    );

    phiHbyA += phig;

    if (adjustFringe)
    {
        oversetAdjustPhi(phiHbyA, U, zoneIdMass);
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U =
                cellMask*
                (
                   HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf)
                );

            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
        if (oversetPatchErrOutput)
        {
            oversetPatchPhiErr(p_rghEqn, phiHbyA);
        }
    }

    #include "continuityErrs.H"

    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf += n*(phi/mesh.magSf() - (n & Uf));
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    // Zero faces H-I for transport Eq after pEq
    phi *= faceMask;

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p_rgh += dimensionedScalar
        (
            "p_rgh",
            p.dimensions(),
            pRefValue - getRefCellValue(p_rgh, pRefCell)
        );
        p == p_rgh + rho*gh;
    }
}
