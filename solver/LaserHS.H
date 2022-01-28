{

    dimensionedScalar pi = constant::mathematical::pi;
    dimensionedScalar a_cond("a_cond",dimensionSet(0, 1, 0, 0, 0),HS_a);//
    dimensionedScalar b_g("b_g",dimensionSet(0, 1, 0, 0, 0),HS_bg);//
    dimensionedScalar v_arc("v_arc",dimensionSet(0, 1, -1, 0, 0),HS_velocity);//
    dimensionedScalar Q_cond("Q_cond",dimensionSet(1, 2, -3, 0, 0),HS_Q);//
    dimensionedScalar lg("lg",dimensionSet(0, 1, 0, 0, 0),HS_lg);//
    dimensionedScalar time("time",dimensionSet(0, 0, 1, 0, 0),runTime.time().value());//



//Set the boundary where the laser is incident
Laser_boundary*=0.0;
Laser_boundary=fvc::average(Laser_boundary);
//Set the boundary where the laser is incident


scalar listLength(0);

DynamicList<vector> initial_points(listLength,vector::zero);
initial_points.clear();//clear the labels

 //List with size equal to number of processors
 List< pointField > gatheredData1(Pstream::nProcs());




Deposition*=0.0;
Testval*=0.0;
forAll( mesh.C(), celli)
{
scalar x_coord = mesh.C()[celli].x();
// scalar y_coord = mesh.C()[celli].y();
scalar z_coord = mesh.C()[celli].z();

// Info<<"x "<<x_coord<<", y "<<y_coord<<", z "<<z_coord<<endl;

if(((Foam::pow(x_coord-b_g.value(),2.0))+(Foam::pow(z_coord-(lg.value()+(v_arc.value()*time.value())),2.0))<=Foam::pow(a_cond.value(),2.0))
&&(Laser_boundary[celli]>SMALL)){
Testval[celli]=1.0;
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
int countbeams=0;
   forAll( pointslistGlobal1, i){
Info<<"Beam: "<<i<<endl;
    //    forAll( mesh.C(), celli)
    //     {
    //         Deposition[celli]=0.0;
    //     }

// vector V2(0, 0.70710678118, 0.70710678118); //ray initial vector at boundary
// vector V2(0, 1, 0); //ray initial vector at boundary
vector V2(0, 0.9, 0.1); //ray initial vector at boundary

    //    Info<<i<<"\t"<<pointslistGlobal1[i]<<endl;


       point V1_tip(pointslistGlobal1[i]);
    //    Foam::sqrt(Foam::pow(pointslistGlobal1[i].x()-b_g.value(),2.0)+Foam::pow(pointslistGlobal1[i].z()-(lg.value()+(v_arc.value()*time.value())),2.0));
    //   scalar r= Foam::sqrt(Foam::pow(pointslistGlobal1[i].x()-b_g.value(),2.0)+Foam::pow(pointslistGlobal1[i].z()-(lg.value()+(v_arc.value()*time.value())),2.0));

    //   scalar Q=((3.0*Q_cond.value())/(a_cond.value()*a_cond.value()*pi.value()))
	// 		*Foam::exp(-3.0*(Foam::pow(((pointslistGlobal1[i].x()-b_g.value())/(a_cond.value())),2.0)+
    //         Foam::pow((pointslistGlobal1[i].z()-(v_arc.value()*time.value())-lg.value())/(a_cond.value()),2.0)));

        scalar Q=1.0;//
       do{
           
            
            // Pout<<V1_tip<<endl;
            point DUMMYMAX(-GREAT,-GREAT,-GREAT);
            scalar DUMMYSCAL(-GREAT);
            // label DUMMYMAXL(-GREAT);

            // List< labelField > Listofcelllabels(Pstream::nProcs());
            // List< pointField > Listofpoints(Pstream::nProcs());

            // List<label> nextlabel(1,-1);
            // List<vector> nextpoint(1,vector::zero);

            label myCellId=mesh.findCell(V1_tip);

            // nextlabel[0]=(mesh.findCell(V1_tip));
            // nextpoint[0]=DUMMYMAX;
            // point V2(-GREAT,-GREAT,-GREAT);
            
                if(myCellId!=-1){

                    
                    Testval[myCellId]=i+1;//set test field to beam flavour

                    if(mag(n_alpha[myCellId])>0.5&&alpha1[myCellId]>0.5){
                            // Q-=1.0;
                        Deposition[myCellId]+=0.5*Q;
                        Q*=0.5;

                            V2=V2-(((((2.0*V2) & n_alpha[myCellId])/(mag(n_alpha[myCellId])*mag(n_alpha[myCellId]))))*n_alpha[myCellId]);//;

                    } else{}

                //////////////////////////////////////////////////////////////Catch rays that get through--maybe gump all their energy here
                if(alpha1[myCellId]>0.5&&mag(n_alpha[myCellId])<0.5){
                    // break;
                    ErrorTrack[myCellId]=1.0;//
                    V2=DUMMYMAX;Q=DUMMYSCAL;
                }//////////////////////////////////////////////////////////////Catch rays that get through


                }
                else{V2=DUMMYMAX;Q=DUMMYSCAL;}

                reduce(V2, maxOp<vector>());//reduce vector //
                reduce(Q, maxOp<scalar>());//
                // reduce(V1_tip, maxOp<vector>());//reduce vector //All processors know where next tip is
                // Pout<<V1_tip<<endl;
                // Pout<<Q<<endl;

                
                // V1_tip+=(iterator_distance*V2);

                label myCellIdnext=mesh.findCell(V1_tip);
                if(myCellIdnext!=-1){
                do{
                    // Pout<<"cell index prior "<<myCellIdnext<<endl;
                    V1_tip+=(iterator_distance*V2);
                     myCellIdnext=mesh.findCell(V1_tip);
                    //  Pout<<"cell index after iteration "<<myCellIdnext<<endl;
                    }
                while(myCellIdnext==myCellId);
                // Pout<<"leaving loop"<<endl;
                    }else{V1_tip=DUMMYMAX;}
                    reduce(V1_tip, maxOp<vector>());//reduce vector //

//  Q-=0.1;
       }while(Q>0);





// countbeams++;
if (countbeams>=1){break;}
   }




// sourceTerm[cellpointhit]=(1.0/yDim[cellpointhit])*










}