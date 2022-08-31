{


    dimensionedScalar time("time",dimensionSet(0, 0, 1, 0, 0),runTime.time().value());//


    double plasma_frequency = Foam::sqrt((e_num_density*constant::electromagnetic::e.value()*constant::electromagnetic::e.value())/(constant::atomic::me.value()*constant::electromagnetic::epsilon0.value()));
    
    double angular_frequency = 2.0*pi.value()*constant::universal::c.value()/wavelength;

    double damping_frequency = plasma_frequency*plasma_frequency*constant::electromagnetic::epsilon0.value()*elec_resistivity;


    double e_r = 1.0-(sqr(plasma_frequency)/(sqr(angular_frequency)+sqr(damping_frequency)));

    double e_i = (damping_frequency/angular_frequency)*((plasma_frequency*plasma_frequency)/((angular_frequency*angular_frequency)+(damping_frequency*damping_frequency)));

    double ref_index = Foam::sqrt((Foam::sqrt((e_r*e_r)+(e_i*e_i))+e_r)/2.0);

    double ext_coefficient = Foam::sqrt((Foam::sqrt((e_r*e_r)+(e_i*e_i))-e_r)/2.0);



    // Info<<"\n"<<endl;
    // Info<<" plasma_frequency \t"<<plasma_frequency<<endl;
    // Info<<" angular_frequency \t"<<angular_frequency<<endl;
    // Info<<" damping_frequency \t"<<damping_frequency<<endl;
    // Info<<" e_r \t"<<e_r<<endl;
    // Info<<" e_i \t"<<e_i<<endl;
    // Info<<" ref_index \t"<<ref_index<<endl;
    // Info<<" ext_coefficient \t"<<ext_coefficient<<endl;
    // Info<<"\n"<<endl;


// vector TESTV1(0, 1, 0); //ray initial vector at boundary
// vector TESTV2(0, 1, 0); //ray initial vector at boundary
// Info<<"\n"<<endl;
// Info<<"testvec"<<Foam::acos((TESTV1 & TESTV2)/(mag(TESTV1)*mag(TESTV2)))<<endl;
// Info<<"\n"<<endl;

//Set the boundary where the laser is incident
Laser_boundary*=0.0;
Laser_boundary=fvc::average(Laser_boundary);
//Set the boundary where the laser is incident


vector normal_interface(0,1,0);//laser comes in on top y boundary

double beam_radius = a_cond.value()/Foam::cos(Foam::acos((normal_interface & (V_incident/mag(V_incident)))/(mag(normal_interface)*mag(V_incident/mag(V_incident)))));
// Adjust sample radius for if beam is not normal too top boundary
double CosTheta_incident = Foam::cos(Foam::acos((normal_interface & (V_incident/mag(V_incident)))/(mag(normal_interface)*mag(V_incident/mag(V_incident)))));
Info<<"cos (theta): \t"<<CosTheta_incident<<endl;

scalar listLength(0);

DynamicList<vector> initial_points(listLength,vector::zero);
initial_points.clear();//clear the labels

 //List with size equal to number of processors
 List< pointField > gatheredData1(Pstream::nProcs());


ErrorTrack*=0.0;
//ErrorTrack2*=0.0;

Deposition*=0.0;
Ray_Number*=0.0;
Ray_Q*=0.0;

forAll( mesh.C(), celli)
{
scalar x_coord = mesh.C()[celli].x();
// scalar y_coord = mesh.C()[celli].y();
scalar z_coord = mesh.C()[celli].z();

// Info<<"x "<<x_coord<<", y "<<y_coord<<", z "<<z_coord<<endl;

// double beam_radius_adjusted_for_initial_incidence_angle = a_cond.value()/Foam::cos(Foam::acos(((V_incident/mag(V_incident)) & normal_interface)/(mag(V_incident/mag(V_incident))*mag(normal_interface))));

if(((Foam::pow(x_coord-b_g.value(),2.0))+(Foam::pow(z_coord-(lg.value()+(v_arc.value()*time.value())),2.0))<=Foam::pow(1.5*beam_radius,2.0))
&&(Laser_boundary[celli]>SMALL)){
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



scalar iterator_distance = (0.5/pi.value())*gMin(yDim);//gMin(xcoord);
Info<<"iterator_distance \t"<<iterator_distance<<endl;
//then loop over all starting points
// int countbeams=0;
   forAll( pointslistGlobal1, i){
Info<<"Beam: "<<i<<endl;

vector V2(V_incident/mag(V_incident));
// Pout<<"TEST_HERE2"<<endl;
    //    Info<<i<<"\t"<<pointslistGlobal1[i]<<endl;


       point V1_tip(pointslistGlobal1[i]);
    
    point mid(b_g.value(),pointslistGlobal1[i].y(),lg.value()+(v_arc.value()*time.value()));
    vector x1 = mid-(10.0*V2);
    vector x2 = mid+(10.0*V2);

    vector x0(pointslistGlobal1[i].x(),pointslistGlobal1[i].y(),pointslistGlobal1[i].z());

    double dist=mag(((x0-x1) ^ (x0-x2)))/mag(x2-x1);//cross product to find distance to beam central axis

    //   scalar Q=((3.0*Q_cond.value())/(a_cond.value()*a_cond.value()*pi.value()))
	// 		*Foam::exp(-3.0*(Foam::pow(((pointslistGlobal1[i].x()-b_g.value())/(beam_radius)),2.0)+
    //         Foam::pow((pointslistGlobal1[i].z()-(v_arc.value()*time.value())-lg.value())/(beam_radius),2.0)));

    scalar Q = CosTheta_incident*((3.0*Q_cond.value())/(Foam::pow(a_cond.value(),2.0)*pi.value()))*Foam::exp(-3.0*((Foam::pow(dist,2.0))/(Foam::pow(a_cond.value(),2.0))));

        // scalar Q=1.0;//
       while(Q>1.0e-9){
           
            
            // Pout<<V1_tip<<endl;
            point DUMMYMAX(-GREAT,-GREAT,-GREAT);
            scalar DUMMYSCAL(-GREAT);
            label myCellId=mesh.findCell(V1_tip);

            
                if(myCellId!=-1){

                    
                    Ray_Number[myCellId]=i+1;//set test field to beam flavour
                    Ray_Q[myCellId]=Q;

                    if(mag(n_filtered[myCellId])>0.5&&alpha_filtered[myCellId]>=0.5){
                        

                            // for(double theta_in=0.0;theta_in<=1.57;theta_in+=0.01){ // to plot absorptivity as a function of incideince angle - a bit hacky
                             double argument = (V2 & n_filtered[myCellId])/(mag(V2)*mag(n_filtered[myCellId]));
                                if(argument>=1.0-SMALL){argument=1.0;}
                                if(argument<=-1.0+SMALL){argument=-1.0;}

                             double theta_in = (std::acos(argument));

                            double alpha_laser = Foam::sqrt((Foam::sqrt(sqr(sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))+(4.0*sqr(ref_index)*sqr(ext_coefficient)))+sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))/(2.0));
                            double beta_laser = Foam::sqrt((Foam::sqrt(sqr(sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))+(4.0*sqr(ref_index)*sqr(ext_coefficient)))-sqr(ref_index)+sqr(ext_coefficient)+sqr(Foam::sin(theta_in)))/(2.0));      
                            double R_s = ((sqr(alpha_laser)+sqr(beta_laser)-(2.0*alpha_laser*Foam::cos(theta_in))+sqr(Foam::cos(theta_in)))/(sqr(alpha_laser)+sqr(beta_laser)+(2.0*alpha_laser*Foam::cos(theta_in))+sqr(Foam::cos(theta_in))));
                            double R_p = R_s*((sqr(alpha_laser)+sqr(beta_laser)-(2.0*alpha_laser*Foam::sin(theta_in)*Foam::tan(theta_in))+(sqr(Foam::sin(theta_in))*sqr(Foam::tan(theta_in))))/(sqr(alpha_laser)+sqr(beta_laser)+(2.0*alpha_laser*Foam::sin(theta_in)*Foam::tan(theta_in))+(sqr(Foam::sin(theta_in))*sqr(Foam::tan(theta_in)))));
                            double absorptivity = 1.0-((R_s+R_p)/2.0);//1.0;//


                                // }

                            if(theta_in>=(pi.value()/2.0)/*-SMALL*/){//sometimes the ray can be reflected and 'skip' along the interface cells - this is unphysical and the ray should traverse  without depositing any energy so set Q to 0 in this instance
                              Q*=0.0;//
                            // Pout<<"vector: "<<V2<<"\t normal: "<<n_filtered[myCellId]<<endl;
                            Deposition[myCellId]+=(absorptivity*Q)/yDim[myCellId];
                            ErrorTrack[myCellId]-=1.0;//
                            // ErrorTrack2[myCellId]=1.0;//
                            //Pout<<"theta \t"<<theta_in<<", absorptivity \t"<<absorptivity<<endl;
                            // V2=-V2;//DUMMYMAX;
                            }
                            // else{}
                            else{
                            // Pout<<"TEST_HERE_pre_reflec"<<endl;
                            Deposition[myCellId]+=(absorptivity*Q)/yDim[myCellId];
                            Q*=(1.0-absorptivity);
                            V2=V2-(((((2.0*V2) & n_filtered[myCellId])/(mag(n_filtered[myCellId])*mag(n_filtered[myCellId]))))*n_filtered[myCellId]);//;
                            // Pout<<"TEST_HERE_reflec"<<endl;
                            }
                    } else{

                        if(alpha_filtered[myCellId]>0.5&&mag(n_filtered[myCellId])<0.5){// if the ray step size happens to be large enough that it skips through the interface send ray back the way it came
                    ErrorTrack[myCellId]+=1.0;//

                            double theta_in = 0.0;//Foam::acos((V2 & n_filtered[myCellId])/(mag(V2)*mag(n_filtered[myCellId])));

                            double alpha_laser = Foam::sqrt((Foam::sqrt(sqr(sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))+(4.0*sqr(ref_index)*sqr(ext_coefficient)))+sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))/(2.0));
                            double beta_laser = Foam::sqrt((Foam::sqrt(sqr(sqr(ref_index)-sqr(ext_coefficient)-sqr(Foam::sin(theta_in)))+(4.0*sqr(ref_index)*sqr(ext_coefficient)))-sqr(ref_index)+sqr(ext_coefficient)+sqr(Foam::sin(theta_in)))/(2.0));    
                            double R_s = ((sqr(alpha_laser)+sqr(beta_laser)-(2.0*alpha_laser*Foam::cos(theta_in))+sqr(Foam::cos(theta_in)))/(sqr(alpha_laser)+sqr(beta_laser)+(2.0*alpha_laser*Foam::cos(theta_in))+sqr(Foam::cos(theta_in))));
                            double R_p = R_s*((sqr(alpha_laser)+sqr(beta_laser)-(2.0*alpha_laser*Foam::sin(theta_in)*Foam::tan(theta_in))+(sqr(Foam::sin(theta_in))*sqr(Foam::tan(theta_in))))/(sqr(alpha_laser)+sqr(beta_laser)+(2.0*alpha_laser*Foam::sin(theta_in)*Foam::tan(theta_in))+(sqr(Foam::sin(theta_in))*sqr(Foam::tan(theta_in)))));


                            double absorptivity = 1.0-((R_s+R_p)/2.0);

                    V2=-V2;//if the ray slips through the interface (unlikely) send it back the way it came because it must have been at 0 degrees anyway
                    //Q=DUMMYSCAL;
                    Deposition[myCellId]+=(absorptivity*Q)/yDim[myCellId];
                        Q*=(1.0-absorptivity);
                } 
                else{}//////////////////////////////////////////////////////////////Catch rays that get through

                    }//

                //////////////////////////////////////////////////////////////Catch rays that get through--maybe gump all their energy here
                


                }
                else{V2=DUMMYMAX;Q=DUMMYSCAL;}

                reduce(V2, maxOp<vector>());//
                reduce(Q, maxOp<scalar>());//


// Pout<<mag(V2)<<endl;
                
                // V1_tip+=(iterator_distance*V2);
                label myCellIdnext=mesh.findCell(V1_tip);
                if(myCellIdnext!=-1){
                while(myCellIdnext==myCellId)
                {
                    //  Pout<<"cell index prior "<<myCellIdnext<<endl;
                    V1_tip+=(iterator_distance*V2);
                     myCellIdnext=mesh.findCell(V1_tip);
                    //   Pout<<"cell index after iteration "<<myCellIdnext<<endl;
                    }
                ;
                // Pout<<"leaving loop"<<endl;
                    }else{V1_tip=DUMMYMAX;}
                    reduce(V1_tip, maxOp<vector>());//reduce vector //

//  Q-=0.1;
       };





// countbeams++;
// if (countbeams>=1){break;}
   }


double TotalQ=gSum(Deposition)*mesh.V()[1];//assumes all cells in mesh are same volume - this code only works on regular meshes anyway
Info<<"Total Q deposited this timestep: \t"<<TotalQ<<endl;
// sourceTerm[cellpointhit]=(1.0/yDim[cellpointhit])*










}