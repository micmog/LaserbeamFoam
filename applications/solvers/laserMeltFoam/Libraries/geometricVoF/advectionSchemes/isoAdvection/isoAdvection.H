/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2019 OpenCFD Ltd.
    Modified code Copyright (C) 2019 DLR
    Modified code Copyright (C) 2018, 2021 Johan Roenby
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::isoAdvection

Description
    An implementation of the isoAdvector geometric Volume-of-Fluid method
    advancing the provided volume fraction field (alpha1) in time using the
    given velocity field, U, and corresponding face fluxes, phi.

    References:

    Main isoAdvector idea:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405
    \endverbatim

    Calculation of rhoPhi:
    \verbatim
        Roenby, J., Bredmose, H., & Jasak, H. (2019).
        IsoAdvector: Geometric VOF on general meshes.
        OpenFOAM® (pp. 281-296). Springer, Cham.
    \endverbatim

    Extension to porous media flows:
    \verbatim
        Missios, K., Jacobsen, N. G., Moeller, K., & Roenby, J.
        Using the isoAdvector Geometric VOF Method for Interfacial Flows
        Through Porous Media. MARINE 2021.
    \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)
    Modified Henning Scheufler, DLR

SourceFiles
    isoAdvection.C
    isoAdvectionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef isoAdvection_H
#define isoAdvection_H

#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "className.H"
#include "reconstructionSchemes.H"
#include "cutFaceAdvect.H"
#include "bitSet.H"
#include "zeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class isoAdvection Declaration
\*---------------------------------------------------------------------------*/

class isoAdvection
{
    // Private data types

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;


    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Dictionary for isoAdvection controls
        const dictionary dict_;

        //- VOF field
        volScalarField& alpha1_;

        //- Often used reference to alpha1 internal field
        scalarField& alpha1In_;

        //- Reference to flux field
        const surfaceScalarField& phi_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Face volumetric water transport
        surfaceScalarField dVf_;

        //- Face volumetric transport
        surfaceScalarField alphaPhi_;

        //- Time spent performing interface advection
        scalar advectionTime_;

        //- store timeIndex
        label timeIndex_;

        // Switches and tolerances. Tolerances need to go into toleranceSwitches

            //- Number of alpha bounding steps
            label nAlphaBounds_;

            //- Tolerance for search of isoFace giving specified VOF value
            scalar isoFaceTol_;

            //- Tolerance for marking of surface cells:
            //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
            scalar surfCellTol_;

            //- Print isofaces in a <case>/isoFaces/isoFaces_#N.vtk files.
            //  Intended for debugging
            bool writeIsoFacesToFile_;

        // Cell and face cutting

            //- List of surface cells
            DynamicLabelList surfCells_;

            //- Pointer to reconstruction scheme
            autoPtr<reconstructionSchemes> surf_;

            //- An isoCutFace object to get access to its face cutting functionality
            cutFaceAdvect advectFace_;

            //- Storage for boundary faces downwind to a surface cell
            DynamicLabelList bsFaces_;

            //- Storage for boundary surface iso face centre
            DynamicVectorList bsx0_;

            //- Storage for boundary surface iso face normal
            DynamicVectorList bsn0_;

            //- Storage for boundary surface iso face speed
            DynamicScalarList bsUn0_;

            //- Switch telling if porosity is enabled
            bool porosityEnabled_;

            //- Pointer to porosity field (read from objectRegistry if exists)
            volScalarField* porosityPtr_;

        // Additional data for parallel runs

            //- List of processor patch labels
            DynamicLabelList procPatchLabels_;

            //- For each patch if it is a processor patch this is a list of the
            //  face labels on this patch that are downwind to a surface cell.
            //  For non-processor patches the list will be empty.
            List<DynamicLabelList> surfaceCellFacesOnProcPatches_;


    // Private Member Functions

        //- No copy construct
        isoAdvection(const isoAdvection&) = delete;

        //- No copy assignment
        void operator=(const isoAdvection&) = delete;


        // Advection functions

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(bitSet& markedCell) const;

            void setProcessorPatches();

            //- For each face calculate volumetric face transport during dt
            void timeIntegratedFlux();

            //- For a given cell return labels of faces fluxing out of this cell
            //  (based on sign of phi)
            void setDownwindFaces
            (
                const label celli,
                DynamicLabelList& downwindFaces
            ) const;

            // LimitFluxes
            template < class SpType, class SuType >
            void limitFluxes
            (
                const SpType& Sp,
                const SuType& Su
            );

            // Bound fluxes
            template < class SpType, class SuType >
            void boundFlux
            (
                const bitSet& nextToInterface,
                surfaceScalarField& dVfcorrectionValues,
                DynamicLabelList& correctedFaces,
                const SpType& Sp,
                const SuType& Su
            );

            //- Given the face volume transport dVf calculates the total volume
            //  leaving a given cell. Note: cannot use dVf member because
            //  netFlux is called also for corrected dVf
            scalar netFlux
            (
                const surfaceScalarField& dVf,
                const label celli
            ) const;

            //- Determine if a cell is a surface cell
            bool isASurfaceCell(const label celli) const
            {
                return
                (
                    surfCellTol_ < alpha1In_[celli]
                 && alpha1In_[celli] < 1 - surfCellTol_
                );
            }

            //- Clear out isoFace data
            void clearIsoFaceData()
            {
                surfCells_.clear();
                bsFaces_.clear();
                bsx0_.clear();
                bsn0_.clear();
                bsUn0_.clear();

            }

        // Face value functions needed for random face access where the face
        // can be either internal or boundary face

            //- Return face value for a given Geometric surface field
            template<typename Type>
            Type faceValue
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei
            ) const;

            //- Set face value for a given Geometric surface field
            template<typename Type>
            void setFaceValue
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei,
                const Type& value
            ) const;


        // Parallel run handling functions

            //- Synchronize dVf across processor boundaries using upwind value
            DynamicList<label> syncProcPatches
            (
                surfaceScalarField& dVf,
                const surfaceScalarField& phi,
                bool returnSyncedFaces=false
            );

            //- Check if the face is on processor patch and append it to the
            //  list of surface cell faces on processor patches
            void checkIfOnProcPatch(const label facei);

            //- Apply the bounding based on user inputs
            void applyBruteForceBounding();

public:

    //- Runtime type information
    TypeName("isoAdvection");

    //- Constructors

        //- Construct given alpha, phi and velocity field. Note: phi should be
        // divergence free up to a sufficient tolerance
        isoAdvection
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~isoAdvection() = default;


    // Member functions

        //- Advect the free surface. Updates alpha field, taking into account
        //  multiple calls within a single time step.
        template < class SpType, class SuType >
        void advect(const SpType& Sp, const SuType& Su);

        // Access functions

            //- Return reconstructionSchemes
            reconstructionSchemes& surf() noexcept
            {
                return surf_();
            }

            //- Return alpha field
            const volScalarField& alpha() const noexcept
            {
                return alpha1_;
            }

            //- Return the controls dictionary
            const dictionary& dict() const noexcept
            {
                return dict_;
            }

            //- Return cellSet of surface cells
            void writeSurfaceCells() const;

            //- Return mass flux
            tmp<surfaceScalarField> getRhoPhi
            (
                const dimensionedScalar& rho1,
                const dimensionedScalar& rho2
            ) const
            {
                return surfaceScalarField::New
                (
                    "rhoPhi",
                    IOobject::NO_REGISTER,
                    (
                        (rho1 - rho2)*dVf_/mesh_.time().deltaT() + rho2*phi_
                    )
                );
            }

            //- Return mass flux
            tmp<surfaceScalarField> getRhoPhi
            (
                const volScalarField& rho1,
                const volScalarField& rho2
            )
            {
                return surfaceScalarField::New
                (
                    "rhoPhi",
                    IOobject::NO_REGISTER,
                    (
                        fvc::interpolate(rho1 - rho2)*alphaPhi_
                      + fvc::interpolate(rho2)*phi_
                    )
                );
            }


// added
            tmp<surfaceScalarField> getRhoPhi
            (
             const surfaceScalarField& rho1,
             const dimensionedScalar rho2
             ) const
            {
                return tmp<surfaceScalarField>
                (
                 new surfaceScalarField
                 (
                  "rhoPhi",
                  (rho1 - rho2)*dVf_/mesh_.time().deltaT() + rho2*phi_
                  )
                 );
            }

            //- reference to alphaPhi
            const surfaceScalarField& alphaPhi() const noexcept
            {
                return alphaPhi_;
            }

            //- time spend in the advection step
            scalar advectionTime() const noexcept
            {
                return advectionTime_;
            }

            //- Write isoface points to .obj file
            void writeIsoFaces
            (
                const DynamicList<List<point>>& isoFacePts
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "isoAdvectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
