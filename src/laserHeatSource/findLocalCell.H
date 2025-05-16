/*---------------------------------------------------------------------------*\
Function
    findLocalCell

Description
    Function that uses a local search to find a cell given its coordinates.
    If the local search fails then default to an N2 search.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef findLocalCell_H
#define findLocalCell_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Find the cell containing the given point and return its cell index. The
    // local search starts from the provided seed cell. If the local search
    // fails then a global search is performed.
    // Searching across processor boundaries is achieved by defaulting to the N2
    // search.
    // The local search will stop after it have performed maxLocalSearch
    // searches; this is to avoid searching the entire domain when loooking for
    // a point that crossed a processor boundary
    label findLocalCell
    (
        const point& p,
        const label seedCellID,
        const polyMesh& mesh,
        const label maxLocalSearch = 100,
        const Switch debug = false
    )
    {
        // If the seed is out-of-range then use N2 search
        if (seedCellID < 0 || seedCellID >= mesh.nCells())
        {
            const label cellID = mesh.findCell(p);
            if (debug)
            {
                Pout<< "    Seed cell out of range: use N2 search. "
                    << "Seed ID = " << seedCellID
                    << ", Cell ID = " << cellID << endl;
            }
            return cellID;
        }

        // Trivial case: check if the point is in the seed cell
        // The default pointInCell decomposes the cell into tets; this is more
        // robust on general polyhedral meshes but is expensive. As we almost
        // exclusively use structured grids, we will use the primitiveMesh
        // implementation based on the face planes, which is similar to
        // polyMesh::cellDecomposition::FACE_PLANES
        // if (mesh.pointInCell(p, seedCellID))
        // if
        // (
        //     mesh.pointInCell
        //     (
        //         p, seedCellID, polyMesh::cellDecomposition::FACE_PLANES
        //     )
        // )
        if (mesh.primitiveMesh::pointInCell(p, seedCellID))
        {
            if (debug)
            {
                Pout<< "    Found point in seed cell. Seed ID = " << seedCellID
                    << ", Cell ID = " << seedCellID
                    << endl;
            }
            return seedCellID;
        }

        // We will perform the search by checking cell-cells
        // We need two lists for storing cell indices: toCheck, and
        // checked

        // Take a reference to the cell-cells
        const labelListList& cellCells = mesh.cellCells();

        // Create toCheck list to store cells to be checked
        // Insert the seed cell-cell neighbours
        SLList<label> toCheck;
        const labelList& seedCellCells = cellCells[seedCellID];
        forAll(seedCellCells, ccI)
        {
            toCheck.append(seedCellCells[ccI]);
        }

        // Create checked list
        labelHashSet checked;

        // Track the number of local searches
        label localSearchI = 0;

        // Perform local search
        while (toCheck.size() > 0)
        {
            // Get next cell to be checked, remove it from the toCheck list and
            // add it to the checked list
            const label cellID = toCheck.removeHead();

            // Check if the point is in this cell
            if (mesh.pointInCell(p, cellID))
	    {
                if (debug)
                {
                    Pout<< "    local search: found after " << localSearchI
                        << " checks. Seed ID = " << seedCellID
                        << ", Cell ID = " << cellID << endl;
                }
                return cellID;
            }

            // Mark the point as checked
            checked.insert(cellID);

            // Insert all cell-cell neighbours of cellID, if they have not
            // already been checked
            const labelList& curCellCells = cellCells[cellID];
            forAll(curCellCells, ccI)
            {
                if (!checked.found(curCellCells[ccI]))
                {
                    toCheck.append(curCellCells[ccI]);
                }
            }

            if (++localSearchI > maxLocalSearch)
            {
                // Break to the N2 search
                if (debug)
                {
                    Pout<< "    Max number of local searches reached: "
                        << maxLocalSearch << endl;
                }
                break;
            }
        }

        // Resort to N2 search
        const label cellID = mesh.findCell(p);
        if (debug)
        {
            Pout<< "    Defaults to N2 search. Seed ID = " << seedCellID
                << ", Cell ID = " << cellID << endl;
        }
        return cellID;
    }

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
