{

    dimensionedScalar pi = constant::mathematical::pi;
    dimensionedScalar a_cond("a_cond",dimensionSet(0, 1, 0, 0, 0),HS_a);//
    dimensionedScalar b_g("b_g",dimensionSet(0, 1, 0, 0, 0),HS_bg);//
    dimensionedScalar v_arc("v_arc",dimensionSet(0, 1, -1, 0, 0),HS_velocity);//
    dimensionedScalar Q_cond("Q_cond",dimensionSet(1, 2, -3, 0, 0),HS_Q);//
    dimensionedScalar lg("lg",dimensionSet(0, 1, 0, 0, 0),HS_lg);//
    dimensionedScalar time("time",dimensionSet(0, 0, 1, 0, 0),runTime.time().value());//
    dimensionedScalar Osc_Amp("Osc_Amp",dimensionSet(0, 1, 0, 0, 0),Oscillation_Amplitude);//
    dimensionedScalar Osc_freq("Osc_freq",dimensionSet(0, 0, -1, 0, 0),Oscillation_Frequency);//


    sourceTerm*=0.0; //SET SOURCETERM TO ZERO
    scalar Y_R(readScalar(meltingDict.lookup("Y_R")));
    scalar focY(readScalar(meltingDict.lookup("focY")));

    scalar foc_shift_vel(readScalar(meltingDict.lookup("foc_shift_vel")));
    scalar Q_ramp_rate(readScalar(meltingDict.lookup("Q_ramp_rate")));
    scalar tshift(readScalar(meltingDict.lookup("tshift")));

    int Rayn=0;
    for(int z_ray=0;z_ray<num_unique_z_coords;z_ray++)
    {
        for(int x_ray=0;x_ray<num_unique_x_coords;x_ray++)
        {
            double miny=1.0e6;
        
            for(int y_ray=0;y_ray<num_unique_y_coords;y_ray++)
            {
            point pnow(xlist[x_ray], ylist[y_ray], zlist[z_ray]);
            label cellpoint = mesh.findCell(pnow);
            //find location of smallest y coordinate
            
                if(alpha1[cellpoint]>HS_deposition_cutoff)
                { 
                //HS_deposition_cutoff 0.99 for conduction mode, 0.01 for keyhole mode
                miny=std::min(miny,ylist[y_ray]);
                }

            }
            point phit(xlist[x_ray], miny, zlist[z_ray]);
            label cellpointhit = mesh.findCell(phit);

            if(yDim[cellpointhit]>1e-12)
                {
                    if(time.value()>tshift)
                {
                double focnow=(focY-(foc_shift_vel*(time.value()-tshift)));
                double Qnow=max((Q_cond.value()-(Q_ramp_rate*(time.value()-tshift))),0.0);

                sourceTerm[cellpointhit]=(1.0/yDim[cellpointhit])*((2.0*Qnow)/
                    ((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((miny-focnow)/Y_R),2.0)))*(a_cond.value()*Foam::sqrt(1.0+Foam::pow(((miny-focnow)/Y_R),2.0)))*pi.value()))
                    *Foam::exp(-2.0*(
                    Foam::pow(((xcoord[cellpointhit]-b_g.value()-(Osc_Amp.value()*Foam::cos(2.0*pi.value()*Osc_freq.value()*time.value())))/((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((miny-focnow)/Y_R),2.0))))),2.0)
                    +Foam::pow((zcoord[cellpointhit]-(v_arc.value()*time.value())-lg.value()/*-(Osc_Amp.value()*Foam::sin(2.0*pi.value()*Osc_freq.value()*time.value()))*/)/((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((miny-focnow)/Y_R),2.0)))),2.0)));
                    }

                    else
                    {
                sourceTerm[cellpointhit]=(1.0/yDim[cellpointhit])*((2.0*Q_cond.value())/
                    ((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((miny-focY)/Y_R),2.0)))*(a_cond.value()*Foam::sqrt(1.0+Foam::pow(((miny-focY)/Y_R),2.0)))*pi.value()))
                    *Foam::exp(-2.0*(
                    Foam::pow(((xcoord[cellpointhit]-b_g.value()-(Osc_Amp.value()*Foam::cos(2.0*pi.value()*Osc_freq.value()*time.value())))/((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((miny-focY)/Y_R),2.0))))),2.0)
                    +Foam::pow((zcoord[cellpointhit]-(v_arc.value()*time.value())-lg.value()/*-(Osc_Amp.value()*Foam::sin(2.0*pi.value()*Osc_freq.value()*time.value()))*/)/((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((miny-focY)/Y_R),2.0)))),2.0)));
                        }
                    }
                    Rayn+=1;
                }//End loop over all rays x
        }//End loop over all rays z


        sourceTerm.correctBoundaryConditions();

                forAll(mesh.C(),celli)
                {
                vector XYZ = mesh.C()[celli];

                       if(time.value()>tshift)
                        {
                        double focnow=(focY-(foc_shift_vel*(time.value()-tshift)));
                        double Qnow=max((Q_cond.value()-(Q_ramp_rate*(time.value()-tshift))),0.0);

                        BeamProfile[celli]=((2.0*Qnow)/
            ((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((XYZ.y()-focnow)/Y_R),2.0)))*(a_cond.value()*Foam::sqrt(1.0+Foam::pow(((XYZ.y()-focnow)/Y_R),2.0)))*pi.value()))
            *Foam::exp(-2.0*(
            Foam::pow(((XYZ.x()-b_g.value()-(Osc_Amp.value()*Foam::cos(2.0*pi.value()*Osc_freq.value()*time.value())))/((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((XYZ.y()-focnow)/Y_R),2.0))))),2.0)
            +Foam::pow((XYZ.z()-(v_arc.value()*time.value())-lg.value()/*-(Osc_Amp.value()*Foam::sin(2.0*pi.value()*Osc_freq.value()*time.value()))*/)/((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((XYZ.y()-focnow)/Y_R),2.0)))),2.0)));
                        }
                        else
                        {
        BeamProfile[celli]=((2.0*Q_cond.value())/
            ((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((XYZ.y()-focY)/Y_R),2.0)))*(a_cond.value()*Foam::sqrt(1.0+Foam::pow(((XYZ.y()-focY)/Y_R),2.0)))*pi.value()))
            *Foam::exp(-2.0*(
            Foam::pow(((XYZ.x()-b_g.value()-(Osc_Amp.value()*Foam::cos(2.0*pi.value()*Osc_freq.value()*time.value())))/((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((XYZ.y()-focY)/Y_R),2.0))))),2.0)
            +Foam::pow((XYZ.z()-(v_arc.value()*time.value())-lg.value()/*-(Osc_Amp.value()*Foam::sin(2.0*pi.value()*Osc_freq.value()*time.value()))*/)/((a_cond.value()*Foam::sqrt(1.0+Foam::pow(((XYZ.y()-focY)/Y_R),2.0)))),2.0)));
                        }
                    }

                }




