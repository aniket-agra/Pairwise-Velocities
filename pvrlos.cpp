#include <mpi.h>
#include <fstream>
#include <cmath>
#include <iostream>

int main(int argc, char* argv[])
{
	double theta,dtheta,theta_max=3.14;
	dtheta=0.1;
	std::ifstream read_file("pvwr");
	std::ofstream write_file("pvlos.dat");
	int size=0;
	double temp=0;
	while(!read_file.eof())
	{
		read_file>>temp>>temp>>temp>>temp;
		size=size+1;
	}
	read_file.close();
	double r[size],v[size],w[size],pvwr[size];
	int rc=0,vc=0,wc=0;
	std::ifstream read_file("pvwr");
	for(int i=0;i<size;i++)
	{
		read_file>>temp;
		if(rc==0)
		{
			r[rc]=temp;
			rc++;
			vc=0;
		}
		else
		if(temp!=r[rc-1])
		{
			r[rc]=temp;
			rc++;
			vc=0;
		}
		read_file>>temp;
		if(vc==0)
		{
			v[vc]=temp;
			vc++;
			wc=0;
		}
		else
		if(temp!=v[vc-1])
		{
			v[vc]=temp;
			vc++;
			wc=0;
		}
		read_file>>temp;
		if(wc==0)
		w[wc]=temp;
		else
		if(temp!=w[wc-1])
		w[wc]=temp;
		wc++;
		read_file>>pvwr[i];
	}

	rc=0;
	while(r[rc]!=NULL)
	{
		theta=0.0;
		while(theta<=theta_max)
		{
			write_file<<r[rc]<<"\t"<<theta<<"\t";
			for(double vlos=-2000.0;vlos<=2000.0;vlos+=100.0)
			{
				if(theta==0.0)
				{
					double v_hr=vlos;
					

				double phi=0.0;
				while(phi<=6.28)
				{
					
	read_file.close();
	MPI::Init(argc,argv);
	MPI::COMM_WORLD.Bcast(&size,1,MPI::INT,0);
	int num_procs=MPI::COMM_WORLD.Get_size();
	int rank=MPI::COMM_WORLD.Get_rank();
	MPI::COMM_WORLD.Reduce(&prob,&final_prob,4040,MPI::INT,MPI::SUM,0);
	MPI::Finalize();
}
