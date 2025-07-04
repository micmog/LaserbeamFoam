/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Initialise the alpha.metal field for a powder bed multiphase simulation
    based on a "locations" file which contains a list of particle locations
    and radii, e.g. as created by the LIGGGHTS DEM code.

Authors
    Gowthaman Parivendhan
    Petar Cosic
    Philip Cardiff

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "objectRegistry.H"
#include "fvMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "volFields.H"
#include <fstream>

using namespace Foam;

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    Info<< "Time = " << runTime.timeName() << endl;

    label count = 0;
    label number = 0;
    scalarField particleR(0);
    scalarField particleX(0);
    scalarField particleY(0);
    scalarField particleZ(0);
    std::string line;
    std::ifstream loc;

    Info<< nl << "Reading the particle coordinates" << endl;

    // Open file and store coordinates in dynamic arrays
    loc.open(runTime.path()/"constant"/"location");
    if (!loc)
    {
        FatalError
            << nl << "Unable to open file 'location'" << nl
            << exit(FatalError);
    }

    while (std::getline(loc, line))
    {
        if (count == 3)
        {
            std::stringstream stream(line);
            stream >> number;
            particleR.setSize(number, 0.0);
            particleX.setSize(number, 0.0);
            particleY.setSize(number, 0.0);
            particleZ.setSize(number, 0.0);
        }
        else if (count > 8)
        {
            std::stringstream stream(line);
            stream
                >> particleX[count - 9]
                >> particleY[count - 9]
                >> particleZ[count - 9]
                >> particleR[count - 9];
        }

        count = count + 1;
    }

    loc.close();

    const IOdictionary bedPlateDict
    (
        IOobject
        (
            "bedPlateDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar xmin = 0;
    scalar xmax = 0;
    scalar ymin = 0;
    scalar ymax = 0;
    scalar zmin = 0;
    scalar zmax = 0;
  
    const bool bedStatus(bedPlateDict.lookupOrDefault<bool>("Bed", false));

    if (bedStatus)
    {
        Info<< "Reading Bed Plate properties" << endl;
        xmin = readScalar(bedPlateDict.lookup("xmin"));
        xmax = readScalar(bedPlateDict.lookup("xmax"));
        ymin = readScalar(bedPlateDict.lookup("ymin"));
        ymax = readScalar(bedPlateDict.lookup("ymax"));
        zmin = readScalar(bedPlateDict.lookup("zmin"));
        zmax = readScalar(bedPlateDict.lookup("zmax"));
    }
    else
    {
        Info<< "Bed plate is switched off" << endl;
    }

    Info<< "Setting field region values" << nl
        << "Number of particles in domain = " << number << endl;

    // Declare Scalar fields alpha for solid,liquid and gas fractions
    Info<< nl << "Reading the alpha.metal field" << endl;
    volScalarField alphaM
    (
        IOobject
        (
            "alpha.metal",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0),
        "zeroGradient"
    );

    // Loop to calculate soild fraction of each mesh cell
    // This currently uses an O(N^2) approach, which is still fine for
    // relatively large cases; if needed, we can change this approach to use
    // a local search
    // Also, the current approach is based purely binary: the cell centre is
    // either inside or outside the particle; we can do better and calculate a
    // volume fraction: TODO!
    const vectorField& CI = mesh.C();
    scalarField& alphaMI = alphaM;
    forAll(CI, cellI)
    {
        const vector& curC = CI[cellI];

        forAll(particleR, particleI)
        {
            const scalar distance =
                Foam::sqrt
                (
                    Foam::pow(curC[vector::X] - particleX[particleI], 2)
                  + Foam::pow(curC[vector::Y] - particleY[particleI], 2)
                  + Foam::pow(curC[vector::Z] - particleZ[particleI], 2)
                );

            if (distance <= particleR[particleI])
            {
                alphaMI[cellI] = 1.0;
                break;
            }

            if (bedStatus)
            {
                if (curC[vector::X] >= xmin && curC[vector::X] <= xmax)
                {
                    if (curC[vector::Y] >= ymin && curC[vector::Y] <= ymax)
                    {
                        if
                        (
                            curC[vector::Z] + SMALL >= zmin
                         && curC[vector::Z] - SMALL <= zmax
                        )
                        {
                            alphaMI[cellI] = 1.0;
                        }
                    }
                }
            }
        }
    }

    alphaM.correctBoundaryConditions();

    Info<< "Writing " << alphaM.name() << " to time = " << runTime.timeName()
        << endl;
    alphaM.write();

    Info<< "\nEnd" << endl;

    return 0;
}

// ************************************************************************* //
