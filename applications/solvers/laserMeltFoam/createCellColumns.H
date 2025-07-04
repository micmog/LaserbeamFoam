//Author: Philip Cardiff, UCD
//Modifications: Gowthaman Parivendhan, UCD

// Get top wall index
const label topWallID = mesh.boundaryMesh().findPatchID(laserPatchName);
if (topWallID == -1)
{
    FatalError
        << "patch not found" << abort(FatalError);
}

// Take a references for convenience
const polyPatch& topWall = mesh.boundaryMesh()[topWallID];
const cellList& cells = mesh.cells();
const faceList& faces = mesh.faces();
const vectorField& points = mesh.points();
const labelList& own = mesh.faceOwner();
const labelList& nei = mesh.faceNeighbour();
    
// Get list of first set of cells
// const unallocLabelList& topWallFaceCells = topWall.faceCells();
const labelList& topWallFaceCells = topWall.faceCells(); // TODO this might not work

// Prepare list of dynamicLists, where each dynamicList will hold an ordered
// column of cell indices, and there will be one of these lists for each
// face on the top patch
PtrList< DynamicList<label> > cellColumns(topWall.size());

// Approximate size for dynamic list intialisation
const label approxSize = 2*mesh.nCells()/topWall.size();

label totalCells=0;

// Populate each column
forAll(cellColumns, colI)
{
    // Initialise dynamic list using a conservative guess for the approximate
    // size
    cellColumns.set(colI, new DynamicList<label>(approxSize));

    // Take a reference to the first list
    DynamicList<label>& curCol = cellColumns[colI];

    // Add first cell
    curCol.append(topWallFaceCells[colI]);

    // Now we will walk from this first cell to the bottom of the column
    // We will walk from face to opposite face in the current cell, then find
    // the next cell and repeat until we get to the opposite boundary
    label curCellID = topWallFaceCells[colI];
    label curFaceID = topWall.start() + colI; // global face index
    bool reachedOtherSide = false;
    do
    {
        // Take a reference to the current cell
        const cell& curCell = cells[curCellID];
        
        // Get the face at the opposite side of the cell
        // Note: this only makes sense for prism type cells
        const label nextFaceID = curCell.opposingFaceLabel(curFaceID, faces);

        if (nextFaceID == -1)
        {
            FatalError
                << "Cannot find the next face!" << abort(FatalError);
        }

        if (mesh.isInternalFace(nextFaceID))
        {
            // Find the cell at the otherside of the nextFace
            label nextCellID = -1;
            if (own[nextFaceID] == curCellID)
            {
                nextCellID = nei[nextFaceID];
            }
            else
            {
                nextCellID = own[nextFaceID];
            }

            // Add nextCellID to the current column
            curCol.append(nextCellID);

            // Update curCellID and curFaceID
            curCellID = nextCellID;
            curFaceID = nextFaceID;
        }
        else
        {
            // We have reached the boundary and are finished!
            reachedOtherSide = true;
        }
    }
    while (!reachedOtherSide);

    // Release extra unused space in the column 
    curCol.shrink();
    totalCells += curCol.size();
}
   
if (totalCells != mesh.nCells())
{
    FatalError
    << "Sum of Cells in list of columns is not equal to the number of cells in mesh." << abort(FatalError);
}
    
forAll(delta,cellI)
{
    const cell& cell = mesh.cells()[cellI];
    labelList pLabels(cell.labels(faces));
    pointField pLocal(pLabels.size(), vector::zero);
        
    forAll(pLabels, pointI)
    {
        pLocal[pointI] = points[pLabels[pointI]];
    }
        
    delta[cellI] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
}
    delta.correctBoundaryConditions();


    // forAll (delta.boundaryField(), patchI) // TODO
    // {
    //     // Forced patch assignment.  HJ, 1/Aug/2010.setFields for reference
    //     delta.boundaryField()[patchI] ==
    //     delta.boundaryField()[patchI].patchInternalField();
    // }
    
    Info<< "min/max(Delta) = "
    << min(delta).value() << ", " << max(delta).value() << endl;

//Sanity check
/*forAll(cellColumns, colI)
{
    DynamicList<label>& curCol = cellColumns[colI];
    Info<<"Column Number: "<<colI<<", Column size: "<<curCol.size()<<endl;
    Info<<"First Element: "<<curCol[0]<<", Last Element: "<<curCol[curCol.size()-1]<<endl;
}

DynamicList<label>& curCol = cellColumns[0];
    
forAll(curCol,colI)
{
    Info<<"Cell ID: "<<curCol[colI]<<endl;
}*/
