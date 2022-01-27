{

sourceTerm*=0.0;
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

			if(epsilon1[cellpoint]<HS_deposition_cutoff)
			{
				miny=std::min(miny,ylist[y_ray]);
			}
		}

		point phit(xlist[x_ray], miny, zlist[z_ray]);
		label cellpointhit = mesh.findCell(phit);
		if(yDim[cellpointhit]>1e-12)
		{
			sourceTerm[cellpointhit]=(cellarea[cellpointhit]/mesh.V()[cellpointhit])*((3.0*Q_cond.value())/(a_cond.value()*c_cond.value()*pi.value()))
			*Foam::exp(-3.0*(Foam::pow(((xcoord[cellpointhit]-b_g.value())/(a_cond.value())),2.0)+Foam::pow((zcoord[cellpointhit]-(v_arc.value()*time.value())-lg.value())/(c_cond.value()),2.0)));
		}
		Rayn+=1;
	}//End loop over all rays x
}//End loop over all rays z

sourceTerm.correctBoundaryConditions();

}



