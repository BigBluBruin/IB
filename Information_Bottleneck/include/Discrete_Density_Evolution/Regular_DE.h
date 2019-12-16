#include "Information_Bottleneck/IB_kernel.h"
#include "Discrete_Density_Evolution/Probability_Combination_Tool.h"
#include "Information_Bottleneck/stats.h"

class Regular_DE
{
public:
    unsigned int dc;
    unsigned int dv;
    double sigma2;
    unsigned int max_iter;
    unsigned int quantization_size;

public: 
    std::vector<std::vector<double>> vari_representation;
    std::vector<std::vector<double>> check_representation;

public:
    Regular_DE(unsigned int Dc, unsigned int Dv, double Sigma2, unsigned int Max_iter,unsigned int Quantization_size);
    void Discrete_Density_Evolution();

};
