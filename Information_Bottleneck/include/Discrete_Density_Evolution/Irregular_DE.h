#include "Information_Bottleneck/IB_kernel.h"
#include "Discrete_Density_Evolution/Probability_Combination_Tool.h"
#include "Information_Bottleneck/stats.h"
#include <string>

class Irregular_DE
{
public:
    std::vector<double> check_edge_dist;
    std::vector<double> vari_edge_dist;
    std::vector<double> eff_check_edge_dist;
    std::vector<double> eff_vari_edge_dist;
    double puncture_rate;
    double sigma2;
    unsigned int max_iter;
    unsigned int quantization_size;
    double stop_threshold;
    double llr_combination_interval;
    unsigned ib_runtime;
    std::string suffix;

public: 
    std::vector<std::vector<double>> vari_representation;
    std::vector<std::vector<double>> check_representation;
    std::vector<std::vector<double>> vari_threshold;
    std::vector<std::vector<double>> check_threshold;

public:
    Irregular_DE(std::vector<double> Check_edge_dist, std::vector<double> Vari_edge_dist, double Sigma2, unsigned int Max_iter,unsigned int Quantization_size,double Stop_threshold,double Llr_combination_interval,unsigned Ib_runtime, std::string Suffix, std::vector<double> Eff_check_edge_dist,std::vector<double> Eff_vari_edge_dist, double Punc_rate);
    int Discrete_Density_Evolution();
    int Discrete_Density_Evolution_punc();

};
