int num_unique_x_coords=0;
int num_unique_z_coords=0;
int num_unique_y_coords=0;

double coord_tolerence=1e-7;

std::vector<double> xlist;
int onetime_1=0;
double x_val=0.0;
double y_val=0.0;
double z_val=0.0;
Info<<"Searching for unique x-coordinates"<<endl;
forAll( mesh.C(), celli)
{
	vector XYZ = mesh.C()[celli];
	if(onetime_1==0)
	{
		y_val=XYZ.y();
		z_val=XYZ.z();
		onetime_1++;
	}
	if(std::abs(XYZ.y()-y_val)<coord_tolerence&&std::abs(XYZ.z()-z_val)<coord_tolerence)
	{
		Info<<XYZ.x()<<"\t"<<XYZ.y()<<"\t"<<XYZ.z()<<endl;
		if(std::find(xlist.begin(),xlist.end(),XYZ.x())==xlist.end())
		{
			xlist.push_back(XYZ.x());
			num_unique_x_coords+=1;
		}
	}
}

std::vector<double> zlist;
onetime_1=0;
x_val=0.0;
y_val=0.0;
z_val=0.0;
Info<<"Searching for unique z-coordinates"<<endl;
forAll( mesh.C(), celli)
{
	vector XYZ = mesh.C()[celli];
	if(onetime_1==0)
	{
		y_val=XYZ.y();
		x_val=XYZ.x();
		onetime_1++;
	}
	
	if(std::abs(XYZ.y()-y_val)<coord_tolerence&&std::abs(XYZ.x()-x_val)<coord_tolerence)
	{
		Info<<XYZ.x()<<"\t"<<XYZ.y()<<"\t"<<XYZ.z()<<endl;
		if(std::find(zlist.begin(),zlist.end(),XYZ.z())==zlist.end())
		{
			zlist.push_back(XYZ.z());
			num_unique_z_coords+=1;
		}
	}
}

std::vector<double> ylist;
onetime_1=0;
x_val=0.0;
y_val=0.0;
z_val=0.0;
Info<<"Searching for unique y-coordinates"<<endl;

forAll( mesh.C(), celli)
{
	vector XYZ = mesh.C()[celli];
	if(onetime_1==0)
	{
		z_val=XYZ.z();
		x_val=XYZ.x();
		onetime_1++;
	}
	
	if(std::abs(XYZ.z()-z_val)<coord_tolerence&&std::abs(XYZ.x()-x_val)<coord_tolerence)
	{
		Info<<XYZ.x()<<"\t"<<XYZ.y()<<"\t"<<XYZ.z()<<endl;
		if(std::find(ylist.begin(),ylist.end(),XYZ.y())==ylist.end())
		{
			ylist.push_back(XYZ.y());
			num_unique_y_coords+=1;
		}
	}
}

double lowest=1.0e6;
double lowest2=1.0e6;

for(int i=0;i<num_unique_y_coords;i++)
{
	if(lowest>ylist[i])
	{
		lowest=ylist[i];
	}
}

for(int i=0;i<num_unique_y_coords;i++)
{
	if(lowest2>ylist[i]&&ylist[i]!=lowest)
	{
		lowest2=ylist[i];
	}
}

Info<<"x: "<<num_unique_x_coords<<" z: "<<num_unique_z_coords<<"y: "<<num_unique_y_coords<<endl;
Info<<"lowest y co-ordinate: "<<lowest<<endl;
Info<<"2nd lowest y co-ordinate: "<<lowest2<<endl;
