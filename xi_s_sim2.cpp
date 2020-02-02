#include <iostream>
#include <fstream>
#include <string>
#include "Math_funcs.hpp"
#include "Three_D_Probability.hpp"

int main()
{
	int size=0,size_corr=0;
	double b=1.0, temp2;
	std::ifstream read_pvwr("p2d.dat"), read_corr("correlation.dat");
	std::string temp;
	std::getline(read_pvwr,temp);
	while(!read_pvwr.eof())
	{
		read_pvwr>>temp2>>temp2>>temp2>>temp2>>temp2;
		//std::cout<<temp<<"\n";
		size++;
	}
	size--;
	read_pvwr.clear();
	read_pvwr.seekg(std::ios::beg);
	double **pvwr;
	pvwr=new double*[size];
	for(int count=0;count<size;count++)
	pvwr[count]=new double[4];
	std::getline(read_pvwr,temp);
	for(int i=0; i<size-1;i++)
	{
	read_pvwr>>pvwr[i][0]>>pvwr[i][1]>>pvwr[i][2]>>pvwr[i][3]>>temp2;
	}
	read_pvwr.close();
	std::cout<<" First line of pvwr = "<<pvwr[0][0]<<" "<<pvwr[0][1]<<" "<<pvwr[0][2]<<" "<<pvwr[0][3]<<"\n";
	while(!read_corr.eof())
	{
		std::getline(read_corr,temp);
		size_corr++;
	}
	size_corr--;
	read_corr.clear();
	read_corr.seekg(std::ios::beg);
	double corr[size_corr], rsample[size_corr];
	for(int i=0; i<size_corr;i++)
	read_corr>>rsample[i]>>corr[i];
	read_corr.close();
	int i=0;
	while(i<size)
	{
		double sep=pvwr[i][0], sum=0.0;
		int j=i;
		while(j<size && pvwr[j][0]==sep)
		{	
			sum+=pvwr[j][3];
			j++;
		}
		for(int k=i;k<j;k++)
		if(sum!=0.0)
		pvwr[k][3]/=(sum*2*pi);
		i=j;
	}
	for(int i=0;i<size;i++)
	if(isnan(pvwr[i][3]))
	{
	std::cout<<"i="<<i<<" r="<<pvwr[i][0]<<" V="<<pvwr[i][1]<<" W="<<pvwr[i][2]<<"\n";
	break;
	}
	std::ofstream write_file("xi_s_sim2.dat");
	std::cout<<size<<"\n";
	for(double s=4.0;s<=202.0;s+=2.0)
	{
		for(double mu=-1.0;mu<=1.0;mu+=0.1)
		{
			double y_max=350.0, y=-y_max, dy=3, spar=s*mu, sper=s*sqrt(1-mu*mu), sum=0.0;
			while(y<=y_max)
			{
				double r=sqrt(sper*sper+y*y);
				if(fabs(y)>=0.01)
				{
					for(double W=0.0;W<=50.0;W+=1.0)
					for(double phi=0.0;phi<=2*pi;phi+=0.1)
					{
						double V=(r*(spar-y)-W*sper*sin(phi))/y, correl=interp(rsample,corr,r,size_corr);
						int V_idx=(int)floor(V+51.0), W_idx=(int)floor(W), r_idx=(int)floor(r/5.0), row=102*r_idx+51*V_idx+W_idx; //102=length of V, 51=length of W
						if(!(r_idx<0 || V_idx<0 || W_idx<0 || row>size-1))
						sum+=(1+b*b*correl)*pvwr[row][3]*0.1*dy;
					}
				}
				y+=dy;
			}
			
			write_file<<s<<"\t"<<mu<<"\t"<<sum-1<<"\n";
		}
		std::cout<<"Result for s="<<s<<" written!\n";
	}
	for(int i=0;i<size;i++)
	delete[] pvwr[i];
	delete[] pvwr;
	write_file.close();
}
