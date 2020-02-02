#include <iostream>
#include <fstream>
#include <string>
#include "Math_funcs.hpp"
#include "Three_D_Probability.hpp"

int main()
{
	double vlos,dvlos=100.0,vlos_max=2000.0,R=2.0,dR=2.0,R_max=300.0,mu=-1.0,dmu=0.1,dW=100.0,W_max=2000.0;
	std::ifstream read_pvwr("pvwr.dat");
	int size=0,l_v=41,l_w=21;
	std::string temp;
	while(!read_pvwr.eof())
	{
		std::getline(read_pvwr,temp);
		size++;
	}
	read_pvwr.clear();
	read_pvwr.seekg(std::ios::beg);
	size--;
	double pvwr[size][4];
	for(int i=0;i<size;i++)
	read_pvwr>>pvwr[i][0]>>pvwr[i][1]>>pvwr[i][2]>>pvwr[i][3];
	std::cout<<"Size="<<size<<"\n";
	std::cout<<"Length of V="<<l_v<<" and of W="<<l_w<<"\n";
	std::cout<<pvwr[12][0]<<" "<<pvwr[12][1]<<" "<<pvwr[12][2]<<" "<<pvwr[12][3]<<"\n";
	read_pvwr.close();
	std::ofstream write_pvlos("pvlos_za.dat");
	while(R<=R_max)
	{
		int r_ind=(int)floor((R-1.0)/2.0),v_ind,w_ind;
		mu=-1.0;
		while(mu<=1.0)
		{	
			write_pvlos<<R<<"\t"<<mu<<"\t";
			double V,W_dummy,W,p;
			for(vlos=-vlos_max;vlos<=vlos_max;vlos+=dvlos)	
			{
				if(fabs(mu)>=0.1)
				W=0.0;
				else
				W=-W_max;
				p=0.0; 
				while(W<=W_max)
				{
					if(fabs(mu)>=0.1)
					{
						W_dummy=W;
						V=(vlos-W_dummy*sqrt(1-mu*mu))/mu;
					}
					else
					{
						V=W;
						W_dummy=vlos;	
					}
					v_ind=(int)floor((V+2150.0)/100.0);
					w_ind=(int)floor((W_dummy+50.0)/100.0);
					if(v_ind<0 || v_ind>=l_v || w_ind<0 || w_ind>=l_w || r_ind*l_v*l_w>=size)
					p+=0.0;
					else
					p+=100.0*pvwr[r_ind*l_v*l_w+v_ind*l_w+w_ind][3];
					W=W+dW;
				}
				write_pvlos<<p<<"\t";
			}
			write_pvlos<<"\n";
			mu=mu+dmu;
		}
		R=R+dR;
	}
	write_pvlos.close();
}