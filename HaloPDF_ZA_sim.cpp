#include <mpi.h>
#include <fstream>
#include <cmath>
#include <iostream>

int main(int argc, char* argv[])
{
	int size=0;
	std::ifstream read_file("cat1");
	double x[7];
	while(!read_file.eof())
	{
		read_file>>x[0]>>x[1]>>x[2]>>x[3]>>x[4]>>x[5]>>x[6];
		size++;
	}
	read_file.close();
	std::cout<<"Size="<<size<<"\n";
	int prob[4040];
	MPI::Init(argc,argv);
	MPI::COMM_WORLD.Bcast(&size,1,MPI::INT,0);
	int num_procs=MPI::COMM_WORLD.Get_size();
	int rank=MPI::COMM_WORLD.Get_rank();
	for(int j=rank; j<size;j+=num_procs)
	{	
		if(j==rank)
		{
			for(int k=0;k<4040;k++)
			prob[k]=0;	
		}
		double m1,x1,y1,z1,vx1,vy1,vz1;
		int i=1;
		std::ifstream read_file("cat1");
		bool end=read_file.eof();
		while(!end && i<=j)
		{
			read_file>>m1>>x1>>y1>>z1>>vx1>>vy1>>vz1;
			end=read_file.eof();
			i++;
		}
		if(!end)
		{	
			
			read_file>>m1>>x1>>y1>>z1>>vx1>>vy1>>vz1;
			while(!read_file.eof())
			{
				double m2,x2,y2,z2,vx2,vy2,vz2;
				read_file>>m2>>x2>>y2>>z2>>vx2>>vy2>>vz2;
				//code to assign minimum of delta x, delta x+2400 and delta x-2400 to delta x, similarly for y and z
				double del_x=(fabs(x2-x1)<=fabs(x2-x1+2400))*(x2-x1)+(fabs(x2-x1)>fabs(x2-x1+2400))*(x2-x1+2400);
				del_x=(fabs(del_x)<=fabs(x2-x1-2400))*del_x+(fabs(del_x)>fabs(x2-x1-2400))*(x2-x1-2400);
				double del_y=(fabs(y2-y1)<=fabs(y2-y1+2400))*(y2-y1)+(fabs(y2-y1)>fabs(y2-y1+2400))*(y2-y1+2400);
				del_y=(fabs(del_y)<=fabs(y2-y1-2400))*del_y+(fabs(del_y)>fabs(y2-y1-2400))*(y2-y1-2400);
				double del_z=(fabs(z2-z1)<=fabs(z2-z1+2400))*(z2-z1)+(fabs(z2-z1)>fabs(z2-z1+2400))*(z2-z1+2400);
				del_z=(fabs(del_z)<=fabs(z2-z1-2400))*del_z+(fabs(del_z)>fabs(z2-z1-2400))*(z2-z1-2400);			
				//calculating separation
				double sep=sqrt(del_x*del_x+del_y*del_y+del_z*del_z);
				//calculating parallel velocity
				double vpar=((vx2-vx1)*del_x+(vy2-vy1)*del_y+(vz2-vz1)*del_z)/sep;
				int m=floor(sep/(double)5);
				int n=floor((vpar+5050)/(double)100);
				int ind=m*101+n;
				if(ind<4040)
				prob[ind]=prob[ind]+1;
			}
		}
			read_file.close();
	}
	int final_prob[4040];
	MPI::COMM_WORLD.Reduce(&prob,&final_prob,4040,MPI::INT,MPI::SUM,0);
	if(rank==0)
	{
		std::ofstream write_file("catalog_1.dat");
		for(int k=0;k<4040;k++)
		write_file<<2.5*(2*floor(k/101)+1)<<"\t"<<-5000+100*(k%101)<<"\t"<<final_prob[k]<<"\n";
		write_file.close();
	}
	MPI::Finalize();
}
