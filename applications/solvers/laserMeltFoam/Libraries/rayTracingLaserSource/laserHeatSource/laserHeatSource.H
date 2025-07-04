/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    laserHeatSource

Description
    Calculate the deposition heat source field from a laser using ray-tracing.

    Laser parameters are defined via the provided dictionary.

    It is assumed that the laser enters the domain at the boundary with an
    outward-facing unit normal of (0 1 0).

    Laser rays are written as VTKs for visualisation.

Author
    Thomas Flint, The University of Manchester.
    Philip Cardiff, University College Dublin.

SourceFiles
    laserHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef laserHeatSource_H
#define laserHeatSource_H

#include "volFields.H"
#include "dictionary.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class laserHeatSource Declaration
\*---------------------------------------------------------------------------*/

class laserHeatSource
:
    public IOdictionary
{
    // Private data

        // Heat source field
        volScalarField deposition_;

        //- Field indicating the boundary faces where the laser enters the
        //  domain
        volScalarField laserBoundary_;

        //- Field for tracking error codes
        volScalarField errorTrack_;

        //- Field for tracking the latest ray to vist each cell
        volScalarField rayNumber_;

        //- Field for tracking the energy deposited in each cell
        volScalarField rayQ_;

        //- y dimension field
        volScalarField yDim_;

        //- Refine flag field
        volScalarField refineFlag_;

        //- Powder sim flag
        //  Needed for cases with spherical particles or spurious currents will
        //  smear out the particle interface un-physically
        const Switch powderSim_;

        //- Time vs laser position interpolation series
        interpolationTable<vector> timeVsLaserPosition_;

        //- Time vs laser power interpolation series
        interpolationTable<scalar> timeVsLaserPower_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        laserHeatSource(const laserHeatSource&);

        //- Disallow default bitwise assignment
        void operator=(const laserHeatSource&);

public:

    //- Runtime type information
    TypeName("laserHeatSource");

    // Constructors

        //- Construct from the mesh
        laserHeatSource
        (
            const fvMesh& mesh
        );

    // Destructor

        virtual ~laserHeatSource()
        {}


    // Member Functions

        //- Const access to the deposition field
        const volScalarField& deposition() const
        {
            return deposition_;
        }

        //- Non-const access to the yDim field
        volScalarField& yDim()
        {
            return yDim_;
        }

        //- Non-const access to the refineFlag field
        volScalarField& refineFlag()
        {
            return refineFlag_;
        }

        //- Const access to the powderSim flag
        const Switch& powderSim() const
        {
            return powderSim_;
        }

        //- Update the deposition field
        //  The filtered volume frame field (alphaFiltered) and the filtered
        //  interface normal field are required to define the interface
        //  The electrical resistivity field is passed as an argument, allowing
        //  a general laser heat source for multi-component and temperature
        //  dependent resistivity solvers
        void updateDeposition
        (
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
