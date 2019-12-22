#include "Information_Bottleneck/IB_kernel.h"
#include "Discrete_Density_Evolution/Probability_Combination_Tool.h"
#include "Information_Bottleneck/stats.h"
#include <string>

class Irregular
{
public:
    std::vector<double> check_edge_dist;
    std::vector<double> vari_edge_dist;
    unsigned int dc;
    unsigned int dv;
    double sigma2;
    unsigned int max_iter;
    unsigned int quantization_size;
    double stop_threshold;

public: 
    std::vector<std::vector<double>> vari_representation;
    std::vector<std::vector<double>> check_representation;
    std::vector<std::vector<double>> vari_threshold;
    std::vector<std::vector<double>> check_threshold;

public:
    Irregular(unsigned int Dc, unsigned int Dv, double Sigma2, unsigned int Max_iter,unsigned int Quantization_size,double Threshold);
    int Discrete_Density_Evolution();

};
