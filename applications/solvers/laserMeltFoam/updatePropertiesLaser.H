// ray tracing
//update the cell size fields
const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();
const vectorField& CI = mesh.C();
scalarField& yDimI = laser.yDim();
forAll(CI, celli)
{
    // vector XYZ = mesh.C()[celli];
    // xcoord[celli]=XYZ.x();
    // zcoord[celli]=XYZ.z();

    const cell& cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll(pLabels, pointi)
    {
        pLocal[pointi] = pp[pLabels[pointi]];
    }

    // xDim[celli] = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
    yDimI[celli] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
    // zDim[celli] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
}

// xDim.correctBoundaryConditions();
laser.yDim().correctBoundaryConditions();
// zDim.correctBoundaryConditions();
//update the cell size fields

// refineflag=fvc::average(alpha1);
metalaverage = fvc::average(alphaM);

laser.refineFlag() *= 0.0;
laser.refineFlag() += metalaverage;

alpha_filtered= alphaM;
alpha_filtered.correctBoundaryConditions();

// LatentHeat=alpha_filtered*LatentHeat1 + (1.0-alpha_filtered)*LatentHeat2;

volVectorField gradepsilon1(fvc::grad(alpha_filtered));

gradepsilon1.correctBoundaryConditions();

n_filtered = fvc::grad(alpha_filtered)/((mag(gradepsilon1) + deltaN));
