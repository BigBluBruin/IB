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
    std::vector<std::vector<double>> joint_xy(2,std::vector<double>(cardi,0));
    joint_xy[0][0]=0.5*(qfunc(-(min-1)/sigma));
    //std::cout<<"------join_x_y first----"<<std::endl;
    if(joint_xy[0][0]==0)
        std::cout<<"zero"<<std::endl;
    joint_xy[1][0]=0.5*(1-qfunc((min+1)/sigma));
    joint_xy[0][cardi-1]=0.5*qfunc((max-1)/sigma);
    joint_xy[1][cardi-1]=0.5*qfunc((max+1)/sigma);
    for(int ii=1;ii<cardi-1;ii++)
    {
        joint_xy[0][ii]=0.5*(qfunc(-(partition_points[ii]-1)/sigma)-qfunc(-(partition_points[ii-1]-1)/sigma));
        joint_xy[1][ii]=0.5*(qfunc((partition_points[ii-1]+1)/sigma)-qfunc((partition_points[ii]+1)/sigma));
    }
    //ave_joinprob(joint_xy);
    return (joint_xy) ;                                                                                                                                                                                                                                                                                                                                                     
}

std::vector<std::vector<double>> gaussian_disretization2 (double min, double max, int cardi, double sigma21, double sigma22)
{
    /*This function realizes uniform quantization of gaussian*/
    double sigma_1=sqrt(sigma21);
    double sigma_2=sqrt(sigma22);
    std::vector<double> partition_points = linspace(min,max,cardi-1);  
    std::vector<std::vector<double>> joint_xy(2,std::vector<double>(cardi,0));
    joint_xy[0][0]=0.5*(qfunc(-(min-1)/sigma_1));
    joint_xy[1][0]=0.5*(1-qfunc((min+1)/sigma_2));
    joint_xy[0][cardi-1]=0.5*qfunc((max-1)/sigma_1);
    joint_xy[1][cardi-1]=0.5*qfunc((max+1)/sigma_2);
    for(int ii=1;ii<cardi-1;ii++)
    {
        joint_xy[0][ii]=0.5*(qfunc(-(partition_points[ii]-1)/sigma_1)-qfunc(-(partition_points[ii-1]-1)/sigma_1));
        joint_xy[1][ii]=0.5*(qfunc((partition_points[ii-1]+1)/sigma_2)-qfunc((partition_points[ii]+1)/sigma_2));
    }
    //ave_joinprob(joint_xy);
    return (joint_xy) ;                                                                                                                                                                                                                                                                                                                                                     
}

std::vector<std::vector<double>> llr_permutation(std::vector<std::vector<double>> &joint_prob, double permutation_factor)
{
    std::vector<double> llr = llr_cal(joint_prob);
    std::vector<std::vector<double>> new_joint_prob(2);
    double tempt;
    for (unsigned ii = 0; ii < llr.size(); ii++)
    {
        if(llr[ii]<0)
        {
            
            tempt=(joint_prob[0][ii]+joint_prob[1][ii])*exp(llr[ii]-permutation_factor)/(1+exp(llr[ii]-permutation_factor));        
            new_joint_prob[0].push_back(tempt);
            new_joint_prob[1].push_back(1-tempt);
            //std::cout<<tempt<<"  "<<std::endl;
        }
        else
        {
            
            new_joint_prob[0].push_back(new_joint_prob[1][llr.size()-1-ii]);
            new_joint_prob[1].push_back(new_joint_prob[0][llr.size()-1-ii]);
        }
        
    }
    ave_joinprob(new_joint_prob);
    return new_joint_prob;
}