dimensionedScalar deltaN = 1e-8/pow(average(mesh.V()),1.0/3.0);

// dimensionedScalar length_unit("length_unit",dimensionSet(0, 1, 0, 0, 0),1.0);

// dimensionedScalar deltaN2 = (1e-8/pow(average(mesh.V()),1.0/3.0))/length_unit;

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

multiphaseMixture mixture(U, phi);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()
);
rho.oldTime();

// Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


volScalarField TSolidus
    (
        IOobject
        (
            "TSolidus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TSolidus",dimensionSet(0, 0, 0, 1, 0),0.0)
    );

          volScalarField TLiquidus
    (
        IOobject
        (
            "TLiquidus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TLiquidus",dimensionSet(0, 0, 0, 1, 0),1.0)
    );

        volScalarField DC
    (
        IOobject
        (
            "DC",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("DC",dimensionSet(1,-3,-1,0,0),1.0e14)
    );

        volScalarField epsilon1
    (
        IOobject
        (
            "epsilon1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionedScalar("epsilon1",dimensionSet(0,0,0,0,0),0.0)
    );
    
    TSolidus=mixture.Tsolidus();
    TLiquidus=mixture.Tliquidus();
    
    epsilon1=max(min((T-TSolidus)/(TLiquidus-TSolidus),scalar(1)),scalar(0));

         volScalarField epsilon1mask
    (
        IOobject
        (
            "epsilon1mask",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
       dimensionedScalar("epsilon1mask",dimensionSet(0,0,0,0,0),0.0)
    );


volScalarField LatentHeat
    (
        IOobject
        (
            "LatentHeat",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("LatentHeat",dimensionSet(0, 2, -2, 0, 0),0.0)
        //(cp*(Temperature-((fluid.Tliquidus()+fluid.Tsolidus())/(scalar(2)))) + delta1*fluid.Latentheat())
    );

        volScalarField LatentHeatGas
    (
        IOobject
        (
            "LatentHeatGas",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("LatentHeatGas",dimensionSet(0, 2, -2, 0, 0),0.0)
        //(cp*(Temperature-((fluid.Tliquidus()+fluid.Tsolidus())/(scalar(2)))) + delta1*fluid.Latentheat())
    );


volScalarField Tvapourisation
    (
        IOobject
        (
            "Tvapourisation",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Tvapourisation",dimensionSet(0, 0, 0, 1, 0),1.0)
        //(cp*(Temperature-((fluid.Tliquidus()+fluid.Tsolidus())/(scalar(2)))) + delta1*fluid.Latentheat())
    );



    volScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0)
    );


volScalarField cp
    (
        IOobject
        (
            "cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("cp",dimensionSet(0, 2, -2, -1, 0),0.0)
        //delta1*fluid.Cp() + delta2*fluid.Cpsolid()
    );


volScalarField TRHS
(
    IOobject
    (
        "TRHS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedScalar("TRHS",dimensionSet(1, -1, -3, 0, 0),scalar(0.0))
);

volScalarField Tcorr
    (
        IOobject
        (
            "Tcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Tcorr",dimensionSet(0, 0, 0, 1, 0),0.0)
    );

volScalarField beta
    (
        IOobject
        (
            "beta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("beta",dimensionSet(0, 0, 0, -1, 0),0.0)
    );

volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rhok",dimensionSet(0,0,0,0,0),1.0)
    );

                volScalarField DilationError
    (
        IOobject
        (
            "DilationError",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("DilationError",dimensionSet(0,0,-1,0,0),0.0)
    );

    

    volVectorField Marangoni
(
    IOobject
    (
        "Marangoni",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Marangoni",dimensionSet(1, -2, -2, 0, 0),vector::zero)
);

volScalarField pVap
(
    IOobject
    (
        "pVap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
); // Recoil Pressure





IOdictionary phaseProperties
(
    IOobject
    (
        "phaseProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


dimensionedScalar p0("p0",dimensionSet(1,-1,-2,0,0,0,0),phaseProperties.lookup("p0"));//atmospheric pressure



volScalarField Qv
(
    IOobject
    (
        "Qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
);// Evaporation Cooling


volScalarField yDim
(
    IOobject
    (
        "yDim",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("yDim", dimensionSet(0, 1, 0, 0, 0), 1.0)
);

const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

forAll( mesh.C(), celli)
{
    // vector XYZ = mesh.C()[celli];
    // xcoord[celli]=XYZ.x();
    // zcoord[celli]=XYZ.z();

    const cell & cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll (pLabels, pointi)
    {
        pLocal[pointi] = pp[pLabels[pointi]];
    }

    // xDim[celli] = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
    yDim[celli] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
    // zDim[celli] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
}

// xDim.correctBoundaryConditions();
yDim.correctBoundaryConditions();
// zDim.correctBoundaryConditions();




volVectorField n_filtered
(
    IOobject
    (
        "n_filtered",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("n_filtered",dimensionSet(0, 0, 0, 0, 0),vector::zero)
);

volScalarField condensate
(
    IOobject
    (
        "condensate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("condensate",dimensionSet(0, 0, 0, 0, 0),0.0)
);

// volScalarField Ray_Number
// (
//     IOobject
//     (
//         "Ray_Number",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("Ray_Number",dimensionSet(0, 0, 0, -0, 0),-1.0)
// );

// volScalarField Deposition
// (
//     IOobject
//     (
//         "Deposition",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("Deposition",dimensionSet(1, -1, -3, -0, 0),-1.0)
// );

// volScalarField Ray_Q
// (
//     IOobject
//     (
//         "Ray_Q",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("Ray_Q",dimensionSet(1, 0, -3, 0, 0),scalar(0.0))
// );

volScalarField Laser_boundary
(
    IOobject
    (
        "Laser_boundary",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Laser_boundary=fvc::average(Laser_boundary);



    // volScalarField e_density
    // (
    //     IOobject
    //     (
    //         "e_density",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("e_density",dimensionSet(0, -3, 0, 0, 0),0.0)
    // );


    volScalarField resistivity
    (
        IOobject
        (
            "resistivity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("resistivity",dimensionSet(1, 3, -3, 0, 0),0.0)
    );

    resistivity=mixture.elec_resistivity();

volScalarField alpha_filtered
(
    IOobject
    (
        "alpha_filtered",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_filtered",dimensionSet(0, 0, 0, -0, 0),0.0)
);


volScalarField molarmass
(
    IOobject
    (
        "molarmass",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
         dimensionedScalar("molarmass", dimensionSet(1, 0, 0, 0, -1), 0.0)//(1, -1, -2, 0, 0)
        // dimensionedScalar("TRHS",dimensionSet(0, 2, -3, 0, 0),scalar(0.0))
);




volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh*rhok
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - rho*gh*rhok;
}

mesh.schemes().setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

// Create laser heat source
laserHeatSource laser(mesh);