/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 DHI
    Copyright (C) 2017 OpenCFD Ltd.
    Copyright (C) 2018 Johan Roenby
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2020 OpenCFD Ltd.
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

Application
    laserMeltFoam

Group
    grpMultiphaseSolvers

Description
    Solver for 2 incompressible, isothermal immiscible fluids and using a VOF
    (volume of fluid) phase-fraction based interface capturing approach and solidification of liquid based on Voller method.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e.  laminar, RAS or LES may be selected.

Author
    Gowthaman Parivendhan, UCD
 
        "When I wrote this, only God and I understood what I was doing
        Now, only God knows"
                                                    -Karl Weierstrass

Ported to OpenFOAM v2412 
     Petar Cosic, UCD 
  
        "Still only God knows"
                                          


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"

// Added 
#include "liquidFractionFunctions.H"
#include "laserThreePhaseMixture.H"
#include "turbulenceModelNew.H"
#include "interpolationTable.H" 
#include "laserHeatSource.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFoam"
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"

    #include "readLaserProperties.H" // added
    #include "createCellColumns.H" //added
    #include "interfaceProperties.H" //added

    #include "initCorrectPhi.H" 
    #include "createUfIfPresent.H"
    #include "porousCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "porousCourantNo.H"
        #include "porousAlphaCourantNo.H"
        #include "setDeltaT.H"
        // #include "createTimeControls.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            threePhaseProperties.correct(); 

            // if (pimple.firstIter() || moveMeshOuterCorrectors) // TODO can we just remove this? mesh will never be moving? 
            // {
            //     // if (isA<dynamicRefineFvMesh>(mesh))
            //     // {
            //     //     advector.surf().reconstruct();
            //     // }

            //     mesh.update();

            //     if (mesh.changing())
            //     {
            //         gh = (g & mesh.C()) - ghRef;
            //         ghf = (g & mesh.Cf()) - ghRef;

            //         // if (isA<dynamicRefineFvMesh>(mesh))
            //         // {
            //         //     advector.surf().mapAlphaField();
            //         //     alpha2 = 1.0 - alpha1;
            //         //     alpha2.correctBoundaryConditions();
            //         //     rho == alpha1*rho1 + alpha2*rho2;
            //         //     rho.correctBoundaryConditions();
            //         //     rho.oldTime() = rho;
            //         //     alpha2.oldTime() = alpha2;
            //         // }

            //         MRF.update();

            //         if (correctPhi)
            //         {
            //             // Calculate absolute flux
            //             // from the mapped surface velocity
            //             phi = mesh.Sf() & Uf();

            //             #include "correctPhi.H"

            //             // Make the flux relative to the mesh motion
            //             fvc::makeRelative(phi, U);

            //             // mixture.correct();
            //         }

            //         if (checkMeshCourantNo)
            //         {
            //             #include "meshCourantNo.H"
            //         }
            //     }
            //   }

            ddtGT.oldTime();

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            if (ray_tracing_on)
            {
              #include "updatePropertiesLaser.H"
              laser.updateDeposition(alpha_filtered, n_filtered, electrical_resistivity);
              Power.value() = laserPowerSeries(runTime.time().timeOutputValue());
            }
            else
            {
              #include "updateLaser.H"
            }

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            #include "TEqn.H"

            #include "continuityErrs.H"

            p = pd + rho * gh;
      
            if (pd.needReference())
            {
              p += dimensionedScalar(
                  "p",
                  p.dimensions(),
                  pRefValue - getRefCellValue(p, pdRefCell));
              pd = p - rho * gh;
            }
      
            #include "updateProperties.H"

            if (pimple.turbCorr())
            {
              turbulence->correct(); 
            }
        }

        volVectorField gradT(fvc::grad(T));

        forAll(gT, cellI)
        {
          if (T[cellI] < Ts.value() && ddtGT.oldTime()[cellI] < 0.0 && gT.oldTime()[cellI] > (0.0 + SMALL))
          {
            scalar t_sol = GREAT;
            if (mag(ddtGT.oldTime()[cellI]) > SMALL)
            {
              t_sol = -gT.oldTime()[cellI] / ddtGT.oldTime()[cellI];
            }
    
            if (t_sol < 0.0)
            {
              t_sol = 0.0;
            }
            else if (t_sol > runTime.deltaT().value())
            {
              t_sol = runTime.deltaT().value();
            }
    
            if (alphaM[cellI] > 0.5 - SMALL)
            {
              solidificationTime[cellI] =
                  runTime.value() - runTime.deltaT().value() + t_sol;
              gradTSol[cellI] = gradT[cellI];
            }
          }
    
          if (T[cellI] > Ts.value()) // Reassign values if there is any remelting
          {
            solidificationTime[cellI] = -1;
            gradTSol[cellI] = vector::zero;
          }
        }


        solidificationTime.correctBoundaryConditions();
        gradTSol.correctBoundaryConditions();
        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

        runTime.printExecutionTime(Info);

        if ( Power.value() == 0.0) // stop the simulation if all cells are solidified
        {

          bool solidified = true;

          if (gMax(alphaL) > 1e-6)
          {
            solidified = false;
          }
    
          if (solidified)
          {
            runTime.writeNow();
            Info << "****All cells solidified****" << endl;
    
            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
            break;
          }
        }


    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
