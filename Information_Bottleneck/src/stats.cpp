#include "Information_Bottleneck/stats.h"

double qfunc(double in)
{
    return 0.5*std::erfc(in/sqrt(2));
}

std::vector<std::vector<double>> gaussian_disretization (double min, double max, int cardi, double sigma2)
{
    /*This function realizes uniform quantization of gaussian*/
    double sigma=sqrt(sigma2);
    std::vector<double> partition_points = linspace(min,max,cardi-1);
    for(const auto& a: partition_points)
        std::cout<<a<<"  ";
    std::cout<<partition_points.size()<<std::endl;    
    std::vector<std::vector<double>> joint_xy(2,std::vector<double>(cardi,0));
    joint_xy[0][0]=0.5*(1-qfunc((min-1)/sigma));
    joint_xy[1][0]=0.5*(1-qfunc((min+1)/sigma));
    joint_xy[0][cardi-1]=0.5*qfunc((max-1)/sigma);
    joint_xy[1][cardi-1]=0.5*qfunc((max+1)/sigma);
    for(int ii=1;ii<cardi-1;ii++)
    {
        joint_xy[0][ii]=0.5*(qfunc((partition_points[ii-1]-1)/sigma)-qfunc((partition_points[ii]-1)/sigma));
        joint_xy[1][ii]=0.5*(qfunc((partition_points[ii-1]+1)/sigma)-qfunc((partition_points[ii]+1)/sigma));
    }
    return joint_xy;                                                                                                                                                                                                                                                                                                                                                     
}