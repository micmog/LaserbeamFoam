zeroField Su;
zeroField Sp;

volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phiCN() + mesh.phi())
  : fvc::div(phiCN())
);

divU *= interpolatedCells*cellMask;
