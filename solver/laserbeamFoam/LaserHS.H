// {
    IOdictionary LaserProperties
    (
        IOobject
        (
            "LaserProperties",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );


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

        // scalar beam_radius_adjusted_for_initial_incidence_angle = a_cond.value()/Foam::cos(Foam::acos(((V_incident/mag(V_incident)) & normal_interface)/(mag(V_incident/mag(V_incident))*mag(normal_interface))));

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
            // Ray_Number[celli]=1.0;
            initial_points.append(mesh.C()[celli]);
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

    // Store the list of cell indices where the ray tips are located; these will
    // be used by the the findLocalSearch function when looking for the new tip
    // cell indices
    labelList rayCellIDs(pointslistGlobal1.size(), -1);

    scalar iterator_distance = (0.5/pi.value())*gMin(yDim);//gMin(xcoord);
    Info<<"iterator_distance \t"<<iterator_distance<<endl;
    //then loop over all starting points
    // int countbeams=0;
    forAll(pointslistGlobal1, i)
    {
        Info<<"Beam: "<<i<<endl;

        vector V2(V_incident/mag(V_incident));
        // Pout<<"TEST_HERE2"<<endl;
        //    Info<<i<<"\t"<<pointslistGlobal1[i]<<endl;


        point V1_tip(pointslistGlobal1[i]);

        point mid(bg_effective.value(),pointslistGlobal1[i].y(),lg_effective.value()+(v_arc.value()*time.value()));
        vector x1 = mid-(10.0*V2);
        vector x2 = mid+(10.0*V2);

        vector x0(pointslistGlobal1[i].x(),pointslistGlobal1[i].y(),pointslistGlobal1[i].z());

        scalar dist=mag(((x0-x1) ^ (x0-x2)))/mag(x2-x1);//cross product to find distance to beam central axis

// Info<<"x0: \t"<<x0<<endl;

// Info<<"dist: \t"<<dist<<endl;


        //   scalar Q=((3.0*Q_cond.value())/(a_cond.value()*a_cond.value()*pi.value()))
        //              *Foam::exp(-3.0*(Foam::pow(((pointslistGlobal1[i].x()-b_g.value())/(beam_radius)),2.0)+
        //         Foam::pow((pointslistGlobal1[i].z()-(v_arc.value()*time.value())-lg.value())/(beam_radius),2.0)));

        scalar Q = CosTheta_incident*((Radius_Flavour*Q_cond.value())/(Foam::pow(a_cond.value(),2.0)*pi.value()))*Foam::exp(-Radius_Flavour*((Foam::pow(dist,2.0))/(Foam::pow(a_cond.value(),2.0))));



        // ID of the processor that contains the beam tip
        label tipProcID = -1;

        while (Q > 1.0e-9)
        {
            // Pout<<V1_tip<<endl;
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


                    // }

                    if(theta_in>=(pi.value()/2.0)/*-SMALL*/) //sometimes the ray can be reflected and 'skip' along the interface cells - this is unphysical and the ray should traverse  without depositing any energy so set Q to 0 in this instance
                    {
                        Q *= 0.0;//
                        // Pout<<"vector: "<<V2<<"\t normal: "<<n_filtered[myCellId]<<endl;
                        Deposition[myCellId]+=(absorptivity*Q)/yDim[myCellId];
                        if (debug)
                        {
                        ErrorTrack[myCellId]-=1.0;//
                        }
                        // ErrorTrack2[myCellId]=1.0;//
                        //Pout<<"theta \t"<<theta_in<<", absorptivity \t"<<absorptivity<<endl;
                        // V2=-V2;//DUMMYMAX;
                    }
                    // else{}
                    else
                    {
                        // Pout<<"TEST_HERE_pre_reflec"<<endl;
                        Deposition[myCellId]+=(absorptivity*Q)/yDim[myCellId];
                        Q *= (1.0-absorptivity);
                        V2=V2-(((((2.0*V2) & n_filtered[myCellId])/(mag(n_filtered[myCellId])*mag(n_filtered[myCellId]))))*n_filtered[myCellId]);//;
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

                        V2=-V2;//if the ray slips through the interface (unlikely) send it back the way it came because it must have been at 0 degrees anyway
                        //Q=DUMMYSCAL;
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


            // // Update seed cells for local search
            // rayCellIDs[i] = myCellIdnext;
        };

        // countbeams++;
        // if (countbeams>=1){break;}
    }


    // assumes all cells in mesh are same volume - this code only works on regular meshes anyway
    const scalar TotalQ = gSum(Deposition)*mesh.V()[0];
    Info<<"Total Q deposited this timestep: \t"<<TotalQ<<endl;
// }
