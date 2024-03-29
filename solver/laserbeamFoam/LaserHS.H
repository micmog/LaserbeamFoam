// {
    // PC: move to createFields.H to avoid IO overhead per time-step
    // IOdictionary LaserProperties
    // (
    //     IOobject
    //     (
    //         "LaserProperties",    // dictionary name
    //         runTime.constant(),     // dict is found in "constant"
    //         mesh,                   // registry for the dict
    //         IOobject::MUST_READ,    // must exist, otherwise failure
    //         IOobject::NO_WRITE      // dict is only read by the solver
    //     )
    // );

    label N_sub_divisions(readLabel(LaserProperties.lookup("N_sub_divisions")));

    const scalar HS_a(readScalar(LaserProperties.lookup("HS_a")));
    const scalar HS_bg(readScalar(LaserProperties.lookup("HS_bg")));
    const scalar HS_velocity(readScalar(LaserProperties.lookup("HS_velocity")));
    const scalar HS_lg(readScalar(LaserProperties.lookup("HS_lg")));
    const scalar HS_Q(readScalar(LaserProperties.lookup("HS_Q")));
    const vector V_incident(LaserProperties.lookup("V_incident"));

    const scalar wavelength(readScalar(LaserProperties.lookup("wavelength")));
    const scalar e_num_density(readScalar(LaserProperties.lookup("e_num_density")));
    const scalar elec_resistivity(readScalar(LaserProperties.lookup("elec_resistivity"))); // this is temperature dependent - will include this in future versions


    const dimensionedScalar pi = constant::mathematical::pi;
    const dimensionedScalar a_cond("a_cond",dimensionSet(0, 1, 0, 0, 0),HS_a);
    const dimensionedScalar b_g("b_g",dimensionSet(0, 1, 0, 0, 0),HS_bg);
    const dimensionedScalar v_arc("v_arc",dimensionSet(0, 1, -1, 0, 0),HS_velocity);
    const dimensionedScalar Q_cond("Q_cond",dimensionSet(1, 2, -3, 0, 0),HS_Q);
    const dimensionedScalar lg("lg",dimensionSet(0, 1, 0, 0, 0),HS_lg);

    const scalar oscAmpX(readScalar(LaserProperties.lookup("HS_oscAmpX")));
    const scalar oscAmpZ(readScalar(LaserProperties.lookup("HS_oscAmpZ")));
    const scalar oscFreqX(readScalar(LaserProperties.lookup("HS_oscFreqX")));
    const scalar oscFreqZ(readScalar(LaserProperties.lookup("HS_oscFreqZ")));


    const dimensionedScalar time("time",dimensionSet(0, 0, 1, 0, 0),runTime.time().value());


    const scalar plasma_frequency = Foam::sqrt((e_num_density*constant::electromagnetic::e.value()*constant::electromagnetic::e.value())/(constant::atomic::me.value()*constant::electromagnetic::epsilon0.value()));
    const scalar angular_frequency = 2.0*pi.value()*constant::universal::c.value()/wavelength;
    const scalar damping_frequency = plasma_frequency*plasma_frequency*constant::electromagnetic::epsilon0.value()*elec_resistivity;
    const scalar e_r = 1.0-(sqr(plasma_frequency)/(sqr(angular_frequency)+sqr(damping_frequency)));
    const scalar e_i = (damping_frequency/angular_frequency)*((plasma_frequency*plasma_frequency)/((angular_frequency*angular_frequency)+(damping_frequency*damping_frequency)));
    const scalar ref_index = Foam::sqrt((Foam::sqrt((e_r*e_r)+(e_i*e_i))+e_r)/2.0);
    const scalar ext_coefficient = Foam::sqrt((Foam::sqrt((e_r*e_r)+(e_i*e_i))-e_r)/2.0);

    //needed for cases with spherical particles or spurious currents will smear out the particle interface un-physically
    const Switch PowderSim(LaserProperties.lookupOrDefault<Switch>("PowderSim", false));

    const scalar dep_cutoff(LaserProperties.lookupOrDefault<scalar>("dep_cutoff", 0.5));

    const scalar Radius_Flavour(LaserProperties.lookupOrDefault<scalar>("Radius_Flavour", 2.0));

    const Switch debug(LaserProperties.lookupOrDefault<Switch>("debug", false));
    const Switch useLocalSearch
    (
        LaserProperties.lookupOrDefault<Switch>("useLocalSearch", true)
    );
    Info<< "useLocalSearch: " << useLocalSearch << endl;
    const label maxLocalSearch
    (
        LaserProperties.lookupOrDefault<label>("maxLocalSearch", 100)
    );

    if (debug)
    {
        Info<<"\n" << nl
            <<" plasma_frequency \t" << plasma_frequency << nl
            <<" angular_frequency \t" << angular_frequency << nl
            <<" damping_frequency \t" << damping_frequency << nl
            <<" e_r \t" << e_r << nl
            <<" e_i \t" << e_i << nl
            <<" ref_index \t" << ref_index << nl
            <<" ext_coefficient \t" << ext_coefficient << nl
            <<"\n" << endl;
    }

    // Set the boundary where the laser is incident
    Laser_boundary *= 0.0;
    Laser_boundary = fvc::average(Laser_boundary);


    const vector normal_interface(0,1,0);//laser comes in on top y boundary

    const scalar beam_radius = a_cond.value()/Foam::cos(Foam::acos((normal_interface & (V_incident/mag(V_incident)))/(mag(normal_interface)*mag(V_incident/mag(V_incident)))));
    // Adjust sample radius for if beam is not normal too top boundary
    const scalar CosTheta_incident = Foam::cos(Foam::acos((normal_interface & (V_incident/mag(V_incident)))/(mag(normal_interface)*mag(V_incident/mag(V_incident)))));
    Info<< "cos (theta): \t" << CosTheta_incident<<endl;

    scalar listLength(0);

    DynamicList<vector> initial_points(listLength,vector::zero);
    initial_points.clear();//clear the labels

    //List with size equal to number of processors
    List< pointField > gatheredData1(Pstream::nProcs());


    ErrorTrack *= 0.0;
    //ErrorTrack2 *= 0.0;

    Deposition *= 0.0;
    Ray_Number *= 0.0;
    Ray_Q *= 0.0;

    dimensionedScalar bg_effective = b_g.value() + oscAmpX*sin(2*pi*oscFreqX*time.value());
    dimensionedScalar lg_effective = lg.value() + oscAmpZ*cos(2*pi*oscFreqZ*time.value());

    forAll(mesh.C(), celli)
    {
        const scalar x_coord = mesh.C()[celli].x();
        // const scalar y_coord = mesh.C()[celli].y();
        const scalar z_coord = mesh.C()[celli].z();

        // scalar beam_radius_adjusted_for_initial_incidence_angle =
        // a_cond.value()/Foam::cos(Foam::acos(((V_incident/mag(V_incident)) & normal_interface)/(mag(V_incident/mag(V_incident))*mag(normal_interface))));

        if
        (
            (
                (Foam::pow(x_coord-bg_effective.value(),2.0))
              + (
                    Foam::pow
                    (
                        z_coord-(lg_effective.value()+(v_arc.value()*time.value())),
                        2.0
                    )
                )
             <= Foam::pow(1.5*beam_radius, 2.0)
            )
        && (Laser_boundary[celli] > SMALL)
        )
        {
            // Ray_Number[celli] = 1.0;
            refineflag2[celli]+=0.5;
            // initial_points.append(mesh.C()[celli]);

            for(label Ray_j=0; Ray_j <N_sub_divisions; Ray_j++){
                for(label Ray_k=0; Ray_k <N_sub_divisions; Ray_k++){
                point p_1(mesh.C()[celli].x()-(yDim[celli]/2.0)+((yDim[celli]/(N_sub_divisions+1))*(Ray_j+1)),mesh.C()[celli].y(),mesh.C()[celli].z()-(yDim[celli]/2.0)+((yDim[celli]/(N_sub_divisions+1))*(Ray_k+1)));
                initial_points.append(p_1);
                }
           }

        }
    }


    //  Populate and gather the list onto the master processor.
    gatheredData1[Pstream::myProcNo()] = initial_points;
    Pstream::gatherList(gatheredData1);

    //  Distibulte the data accross the different processors
    Pstream::scatterList(gatheredData1);

    pointField pointslistGlobal1//list of initial points
    (
        ListListOps::combine<Field<vector> >
        (
            gatheredData1,
            accessOp<Field<vector> >()
        )
    );

    // For each beam, store the starting point and locations at which the rays
    // change direction. Also, store the global ordered index of the ray
    // direction-change points
    PtrList<DynamicList<vector>> beamDirectionChangePoints
    (
        pointslistGlobal1.size()
    );
    PtrList<DynamicList<int>> beamDirectionChangeOrder
    (
        pointslistGlobal1.size()
    );

    // Initialise: beams will likely change direction less than 100 times
    forAll(beamDirectionChangePoints, rayI)
    {
        beamDirectionChangePoints.set
        (
            rayI,
            new DynamicList<vector>(100)
        );
        beamDirectionChangeOrder.set
        (
            rayI,
            new DynamicList<int>(100)
        );

        // Add initial point
        beamDirectionChangePoints[rayI].append(pointslistGlobal1[rayI]);
        beamDirectionChangeOrder[rayI].append(0);
    }

    // Store the list of cell indices where the ray tips are located; these will
    // be used by the the findLocalSearch function when looking for the new tip
    // cell indices
    labelList rayCellIDs(pointslistGlobal1.size(), -1);

    scalar iterator_distance = (0.5/pi.value())*gMin(yDim);//gMin(xcoord);
    Info<<"iterator_distance \t"<<iterator_distance<<endl;
    // Loop over all starting points
    forAll(pointslistGlobal1, i)
    {
        Info<< "Beam: "<< i <<endl;

        vector V2(V_incident/mag(V_incident));
        point V1_tip(pointslistGlobal1[i]);
        const point mid
        (
            bg_effective.value(),
            pointslistGlobal1[i].y(),
            lg_effective.value() + (v_arc.value()*time.value())
        );
        const vector x1 = mid - (10.0*V2);
        const vector x2 = mid + (10.0*V2);
        const vector x0
        (
            pointslistGlobal1[i].x(),
            pointslistGlobal1[i].y(),
            pointslistGlobal1[i].z()
        );

        // Cross product to find distance to beam central axis
        const scalar dist = mag(((x0-x1) ^ (x0-x2)))/mag(x2 - x1);

        // Global index to track the order of the ray direction-changes
        // This is only used for post-processing to write VTKs of the beams
        label directionChangeOrderI = 0;

        // Info<<"x0: \t"<<x0<<endl;

// Info<<"dist: \t"<<dist<<endl;


        //   scalar Q=((3.0*Q_cond.value())/(a_cond.value()*a_cond.value()*pi.value()))
        //              *Foam::exp(-3.0*(Foam::pow(((pointslistGlobal1[i].x()-b_g.value())/(beam_radius)),2.0)+
        //         Foam::pow((pointslistGlobal1[i].z()-(v_arc.value()*time.value())-lg.value())/(beam_radius),2.0)));

        scalar Q = (CosTheta_incident/(N_sub_divisions*N_sub_divisions))*((Radius_Flavour*Q_cond.value())/(Foam::pow(a_cond.value(),2.0)*pi.value()))*Foam::exp(-Radius_Flavour*((Foam::pow(dist,2.0))/(Foam::pow(a_cond.value(),2.0))));



        // ID of the processor that contains the beam tip
        label tipProcID = -1;

        while (Q > 1.0e-9)
        {
            // Track when the tip changes direction for post-processing the rays
            bool beamChangedDirection = false;

            point DUMMYMAX(-GREAT,-GREAT,-GREAT);
            scalar DUMMYSCAL(-GREAT);

            // Search for the cell that contains the local beam tip
            // Only the processor that contained the old tip will perform the
            // search, or all processor will search if the old tip is not on any
            // processor
            label myCellId = -1;
            if (tipProcID == Pstream::myProcNo() || tipProcID == -1)
            {
                if (useLocalSearch)
                {
                    myCellId =
                        findLocalCell(V1_tip, rayCellIDs[i], mesh, maxLocalSearch, debug);
                }
                else
                {
                    myCellId = mesh.findCell(V1_tip);
                }
            }

            // Proc ID where the tip is located
            // If the tip in not on any processor, then this is set to -1
            if (myCellId != -1)
            {
                tipProcID = Pstream::myProcNo();
            }
            else
            {
                tipProcID = -1;
            }
            reduce(tipProcID, maxOp<label>());

            if (myCellId != -1)
            {
                Ray_Number[myCellId]=i+1;//set test field to beam flavour
                Ray_Q[myCellId]=Q;

                if(mag(n_filtered[myCellId])>0.5&&alpha_filtered[myCellId]>=dep_cutoff)
                {


                    // for(scalar theta_in=0.0;theta_in<=1.57;theta_in+=0.01){ // to plot absorptivity as a function of incideince angle - a bit hacky
                    scalar argument = (V2 & n_filtered[myCellId])/(mag(V2)*mag(n_filtered[myCellId]));
                    if(argument>=1.0-SMALL)
                    {
                        argument=1.0;
                    }
                    if(argument<=-1.0+SMALL)
                    {
                        argument=-1.0;
                    }

                    scalar theta_in = (std::acos(argument));

                    scalar alpha_laser = Foam::sqrt((Foam::sqrt(sqr(sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))+(4.0*sqr(ref_index)*sqr(ext_coefficient)))+sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))/(2.0));
                    scalar beta_laser = Foam::sqrt((Foam::sqrt(sqr(sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))+(4.0*sqr(ref_index)*sqr(ext_coefficient)))-sqr(ref_index)+sqr(ext_coefficient)+sqr(Foam::sin(theta_in)))/(2.0));
                    scalar R_s = ((sqr(alpha_laser)+sqr(beta_laser)-(2.0*alpha_laser*Foam::cos(theta_in))+sqr(Foam::cos(theta_in)))/(sqr(alpha_laser)+sqr(beta_laser)+(2.0*alpha_laser*Foam::cos(theta_in))+sqr(Foam::cos(theta_in))));
                    scalar R_p = R_s*((sqr(alpha_laser)+sqr(beta_laser)-(2.0*alpha_laser*Foam::sin(theta_in)*Foam::tan(theta_in))+(sqr(Foam::sin(theta_in))*sqr(Foam::tan(theta_in))))/(sqr(alpha_laser)+sqr(beta_laser)+(2.0*alpha_laser*Foam::sin(theta_in)*Foam::tan(theta_in))+(sqr(Foam::sin(theta_in))*sqr(Foam::tan(theta_in)))));
                    scalar absorptivity = 1.0-((R_s+R_p)/2.0);//1.0;//
                    // scalar absorptivity = 1.0;//1.0;//


                    // }

                    // Sometimes the ray can be reflected and 'skip' along the
                    // interface cells - this is unphysical and the ray should
                    // traverse  without depositing any energy so set Q to 0 in
                    // this instance
                    if (theta_in>=(pi.value()/2.0))
                    {
                        Q *= 0.0;
                        Deposition[myCellId]+=(absorptivity*Q)/yDim[myCellId];
                        if (debug)
                        {
                            ErrorTrack[myCellId] -= 1.0;
                        }
                        beamChangedDirection = true;
                    }
                    // else{}
                    else
                    {
                        // Pout<<"TEST_HERE_pre_reflec"<<endl;
                        Deposition[myCellId]+=(absorptivity*Q)/yDim[myCellId];
                        Q *= (1.0-absorptivity);
                        V2=V2-(((((2.0*V2) & n_filtered[myCellId])/(mag(n_filtered[myCellId])*mag(n_filtered[myCellId]))))*n_filtered[myCellId]);//;
                        beamChangedDirection = true;
                        // Pout<<"TEST_HERE_reflec"<<endl;
                    }
                }
                else
                {

                    if(alpha_filtered[myCellId]>dep_cutoff&&mag(n_filtered[myCellId])<0.5) // if the ray step size happens to be large enough that it skips through the interface send ray back the way it came
                    {
                        if (debug)
                        {
                        ErrorTrack[myCellId]+=1.0;//
                        }
                        scalar theta_in = 0.0;//Foam::acos((V2 & n_filtered[myCellId])/(mag(V2)*mag(n_filtered[myCellId])));

                        scalar alpha_laser = Foam::sqrt((Foam::sqrt(sqr(sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))+(4.0*sqr(ref_index)*sqr(ext_coefficient)))+sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))/(2.0));
                        scalar beta_laser = Foam::sqrt((Foam::sqrt(sqr(sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))+(4.0*sqr(ref_index)*sqr(ext_coefficient)))-sqr(ref_index)+sqr(ext_coefficient)+sqr(Foam::sin(theta_in)))/(2.0));
                        scalar R_s = ((sqr(alpha_laser)+sqr(beta_laser)-(2.0*alpha_laser*Foam::cos(theta_in))+sqr(Foam::cos(theta_in)))/(sqr(alpha_laser)+sqr(beta_laser)+(2.0*alpha_laser*Foam::cos(theta_in))+sqr(Foam::cos(theta_in))));
                        scalar R_p = R_s*((sqr(alpha_laser)+sqr(beta_laser)-(2.0*alpha_laser*Foam::sin(theta_in)*Foam::tan(theta_in))+(sqr(Foam::sin(theta_in))*sqr(Foam::tan(theta_in))))/(sqr(alpha_laser)+sqr(beta_laser)+(2.0*alpha_laser*Foam::sin(theta_in)*Foam::tan(theta_in))+(sqr(Foam::sin(theta_in))*sqr(Foam::tan(theta_in)))));


                        scalar absorptivity = 1.0-((R_s+R_p)/2.0);
                        // scalar absorptivity = 1.0;//

                        V2=-V2;//if the ray slips through the interface (unlikely) send it back the way it came because it must have been at 0 degrees anyway
                        //Q=DUMMYSCAL;
                        beamChangedDirection = true;
                        Deposition[myCellId]+=(absorptivity*Q)/yDim[myCellId];
                        Q *= (1.0-absorptivity);
                    }
                    else {} //////////////////////////////////////////////////////////////Catch rays that get through

                }//

                //////////////////////////////////////////////////////////////Catch rays that get through--maybe gump all their energy here



            }
            else
            {
                // The tip is not on this processor for one of two reasons:
                // 1. the tip left the entire global domain
                // 2. the tip is on another processor
                V2 = DUMMYMAX;
                Q = DUMMYSCAL;
                beamChangedDirection = true;
            }

            reduce(V2, maxOp<vector>());
            reduce(Q, maxOp<scalar>());

            // Update seed cells for local search
            rayCellIDs[i] = myCellId;

            // V1_tip+=(iterator_distance*V2);
            //label myCellIdnext=mesh.findCell(V1_tip);
            if (tipProcID == Pstream::myProcNo())
            {
                label myCellIdnext =
                    findLocalCell(V1_tip, rayCellIDs[i], mesh, maxLocalSearch, debug);

                if (myCellIdnext != -1)
                {
                    while (myCellIdnext == myCellId)
                    {
                        if (beamChangedDirection)
                        {
                            // Write current tip position to array
                            beamDirectionChangePoints[i].append(V1_tip);
                            beamDirectionChangeOrder[i].append(directionChangeOrderI);
                            beamChangedDirection = false;
                        }

                        V1_tip += (iterator_distance*V2);
                        //myCellIdnext = mesh.findCell(V1_tip);
                        myCellIdnext =
                            findLocalCell(V1_tip, rayCellIDs[i], mesh, maxLocalSearch, debug);
                    }
                }
                else
                {
                    V1_tip=DUMMYMAX;
                }

                // Update direction-change ordered index
                directionChangeOrderI++;

                // Update seed cells for local search
                rayCellIDs[i] = myCellIdnext;
            }
            else
            {
                V1_tip=DUMMYMAX;
            }
            reduce(V1_tip, maxOp<vector>());//reduce vector //
            //  Q-=0.1;

            if (rayCellIDs[i] == -1)
            {
                tipProcID = -1;
            }
            reduce(tipProcID, maxOp<label>());

            // Sync direction-change ordered index
            reduce(directionChangeOrderI, maxOp<int>());

            // // Update seed cells for local search
            // rayCellIDs[i] = myCellIdnext;
        };

        // countbeams++;
        // if (countbeams>=1){break;}
    }


    //  gSum(mesh.V()*Deposition);
    // const scalar TotalQ = gSum(Deposition*mesh.V().value());
   const scalar TotalQ = fvc::domainIntegrate(Deposition).value();
    Info<< "Total Q deposited this timestep: \t" << TotalQ <<endl;
// }

// Combine rays across procs
if (runTime.outputTime() && Pstream::parRun())
{
    if (debug)
    {
        Info<< "Parallel syncing beams!" << endl;
    }

    // The ray starting points were added to the beamDirectionChangePoints list
    // on all procs, so we will remove them from all procs apart from the master
    // Note: the beamDirectionChangePoints list is only synced at output times
    // and will only be correct on the master proc which writes them
    if (!Pstream::master())
    {
        forAll(beamDirectionChangePoints, rayI)
        {
            beamDirectionChangePoints[rayI] =
                SubField<vector>
                (
                    beamDirectionChangePoints[rayI],
                    beamDirectionChangePoints[rayI].size() - 1,
                    1
                );

            beamDirectionChangeOrder[rayI] =
                SubField<int>
                (
                    beamDirectionChangeOrder[rayI],
                    beamDirectionChangeOrder[rayI].size() - 1,
                    1
                );
        }
    }

    // Sync beams across procs
    forAll(beamDirectionChangePoints, rayI)
    {
        {
            List<List<vector>> gatheredField(Pstream::nProcs());
            gatheredField[Pstream::myProcNo()] = beamDirectionChangePoints[rayI];
            Pstream::gatherList(gatheredField);

            beamDirectionChangePoints[rayI] =
                ListListOps::combine<List<vector>>
                (
                    gatheredField,
                    accessOp<List<vector>>()
                );
        }

        {
            List<List<int>> gatheredField(Pstream::nProcs());
            gatheredField[Pstream::myProcNo()] = beamDirectionChangeOrder[rayI];
            Pstream::gatherList(gatheredField);

            beamDirectionChangeOrder[rayI] =
                ListListOps::combine<List<int>>
                (
                    gatheredField,
                    accessOp<List<int>>()
                );
        }

        // Re-order the list
        if (Pstream::master())
        {
            SortableList<int> sortedOrder(beamDirectionChangeOrder[rayI]);
            List<vector> unsortedPoints(beamDirectionChangePoints[rayI]);
            List<int> unsortedOrder(beamDirectionChangeOrder[rayI]);
            forAll(sortedOrder, i)
            {
                beamDirectionChangePoints[rayI][i] =
                    unsortedPoints[sortedOrder.indices()[i]];
                beamDirectionChangeOrder[rayI][i] =
                    unsortedOrder[sortedOrder.indices()[i]];
            }
        }
    }
}


// Write rays
if (runTime.outputTime() && Pstream::master())
{
    if (debug)
    {
        forAll(beamDirectionChangePoints, rayI)
        {
            Info<< "ray " << rayI << endl;
            forAll(beamDirectionChangePoints[rayI], i)
            {
                Info<< "    " << beamDirectionChangePoints[rayI][i] << endl;
            }
        }
    }

    // Write rays in VTK format
    // See
    // https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html

    // Create a directory for the VTK files
    fileName vtkDir;
    if (Pstream::parRun())
    {
        vtkDir = runTime.path()/".."/"VTKs";
    }
    else
    {
        vtkDir = runTime.path()/"VTKs";
    }

    mkDir(vtkDir);

    // Create a VTK file
    OFstream rayVtkFile
    (
        vtkDir/"rays_" //+ runTime.timeName() + "_"
      + Foam::name(runTime.timeIndex()) + ".vtk"
    );

    Info<< "Writing rays to " << rayVtkFile.name() << endl;

    // Write header
    rayVtkFile
        << "# vtk DataFile Version 2.0" << nl
        << "Rays" << nl
        << "ASCII" << endl;

    // Count the number of points and calculate the offset for each ray
    label nRayPoints = 0;
    labelList pointIdOffset(beamDirectionChangePoints.size(), 0);
    forAll(beamDirectionChangePoints, rayI)
    {
        nRayPoints += beamDirectionChangePoints[rayI].size();

        if (rayI > 0)
        {
            pointIdOffset[rayI] =
                pointIdOffset[rayI - 1]
              + beamDirectionChangePoints[rayI - 1].size();
        }
    }

    // Write points
    rayVtkFile
        << "DATASET POLYDATA" << nl
        << "POINTS " << nRayPoints << " double" << endl;

    // Add ray points
    forAll(beamDirectionChangePoints, rayI)
    {
        forAll(beamDirectionChangePoints[rayI], i)
        {
            rayVtkFile
                << beamDirectionChangePoints[rayI][i].x() << " "
                << beamDirectionChangePoints[rayI][i].y() << " "
                << beamDirectionChangePoints[rayI][i].z() << endl;
        }
    }

    // Count the number of lines
    label nRayLines = 0;
    forAll(beamDirectionChangePoints, rayI)
    {
        // Note: we must add 1 as the VTK format requires it
        nRayLines += beamDirectionChangePoints[rayI].size() + 1;
    }

    // Write lines
    rayVtkFile
        << "LINES " << beamDirectionChangePoints.size() << " " << nRayLines
        << endl;

    forAll(beamDirectionChangePoints, rayI)
    {
        // Write the number of points in the line
        rayVtkFile
            << beamDirectionChangePoints[rayI].size();

        // Write indices of points
        forAll(beamDirectionChangePoints[rayI], i)
        {
            rayVtkFile
                << " " << pointIdOffset[rayI] + i;
        }

        rayVtkFile
            << endl;
    }
}
