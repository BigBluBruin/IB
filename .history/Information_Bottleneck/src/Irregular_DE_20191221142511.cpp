#include "Discrete_Density_Evolution/Irregular_DE.h"

Irregular_DE::Irregular_DE(std::vector<double> Check_edge_dist={0.0}, std::vector<double> Vari_edge_dist={0.0}, double Sigma2=0.0, unsigned int Max_iter=0,unsigned int Quantization_size=0,double Stop_threshold=0.0)
{
    check_edge_dist=Check_edge_dist;
    vari_edge_dist=Vari_edge_dist;
    sigma2=Sigma2;
    max_iter=Max_iter;
    quantization_size=Quantization_size;
    stop_threshold=Stop_threshold;
}


int Irregular_DE::Discrete_Density_Evolution()
{
    std::vector<std::vector<double>> combined_check_dist, combined_vari_dist;
    std::vector<std::vector<double>> channel_observation = gaussian_disretization(-5, 5, 3000, sigma2);
    IB_kernel channel_IB(channel_observation, quantization_size, 2000);
    channel_IB.smIB();
    std::vector<std::vector<double>> first_input = channel_IB.prob_join_xt;
    std::vector<std::vector<double>> second_input = first_input;
    std::vector<std::vector<double>> combined_output, prob_llr_combined;
}