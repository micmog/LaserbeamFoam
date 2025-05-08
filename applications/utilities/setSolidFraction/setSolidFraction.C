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
    Selects a cell set through a dictionary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "objectRegistry.H"
// #include "foamTime.H"
#include "fvMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "volFields.H"
#include <fstream>

using namespace Foam;
using namespace std;

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"

  Info << "Time = " << runTime.timeName() << endl;

#include "createMesh.H"

  int count = 0, number;
  float *r, *x, *y, *z; // Dynamic arrays to store location and radius of the particles
  std::string line;
  ifstream loc;

  Info << "\nReading coordinates of Particles" << endl;

  // Open file and store coordinates in dynamic arrays
  loc.open(runTime.path() / "constant" / "location");
  if (!loc)
  {
    Info << "\nUnable to open file 'location'\n";
    exit(1);
  }

  while (std::getline(loc, line))
  {
    if (count == 3)
    {
      std::stringstream stream(line);
      stream >> number;
      r = new float[number];
      x = new float[number];
      y = new float[number];
      z = new float[number];
    }

    if (count > 8)
    {
      std::stringstream stream(line);
      stream >> x[count - 9] >> y[count - 9] >> z[count - 9] >> r[count - 9];
    }
    count = count + 1;
  }

  loc.close();

  IOdictionary bedPlateDict(
      IOobject(
          "bedPlateDict",
          runTime.system(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE));

  scalar xmin = 0;
  scalar xmax = 0;
  scalar ymin = 0;
  scalar ymax = 0;
  scalar zmin = 0;
  scalar zmax = 0;
  
  const bool bedStatus(bedPlateDict.lookupOrDefault<bool>("Bed", false));

  if (bedStatus)
  {
    Info << "Reading Bed Plate properties" << endl;
    bedPlateDict.lookup("xmin") >> xmin;
    bedPlateDict.lookup("xmax") >> xmax;
    bedPlateDict.lookup("ymin") >> ymin;
    bedPlateDict.lookup("ymax") >> ymax;
    bedPlateDict.lookup("zmin") >> zmin;
    bedPlateDict.lookup("zmax") >> zmax;
  }
  else
  {
    Info << "Bed plate is switched off" << endl;
  }

  Info << "Setting field region values" << endl;
  Info << "\nNumber of Particles in domain :" << number << endl;

  // Declare Scalar fields alpha for solid,liquid and gas fractions
  volScalarField alphaM(
      IOobject(
          "alpha.metal",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE),
      mesh);

  // Initialise field values to 0.0
  forAll(alphaM.internalField(), cellI)
  {
    alphaM[cellI] = 0.0;
  }

  // Declare scalarFields to get x,y and z coords of Mesh cell centres
  volScalarField px(mesh.C().component(0));
  volScalarField py(mesh.C().component(1));
  volScalarField pz(mesh.C().component(2));

  // Loop to calculate Soild fraction of each mesh cell
  forAll(alphaM.internalField(), cellI)
  {
    for (int j = 0; j < number; j++)
    {
      scalar distance = Foam::sqrt(Foam::pow(px[cellI] - (x[j]), 2) + Foam::pow(py[cellI] - (y[j]), 2) + Foam::pow(pz[cellI] - (z[j]), 2));

      if (distance <= (r[j]))
      {
        alphaM[cellI] = 1.0;
        break;
      }

      if (bedStatus)
      {
        if (px[cellI] >= xmin && px[cellI] <= xmax)
        {
          if (py[cellI] >= ymin && py[cellI] <= ymax)
          {
            if (pz[cellI] + SMALL >= zmin && pz[cellI] - SMALL <= zmax)
            {
              alphaM[cellI] = 1.0;
            }
          }
        }
      }
    }
  }

  // forAll (alphaM.boundaryField(), patchi)
  // {
  // // Forced patch assignment.  HJ, 1/Aug/2010.
  // alphaM.boundaryField()[patchi] ==
  // alphaM.boundaryField()[patchi].patchInternalField();
  // }

  alphaM.write();

  // Always clean up when you are finished
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] r;
  return 0;
  Info << "\nEnd" << endl;
}

// ************************************************************************* //
