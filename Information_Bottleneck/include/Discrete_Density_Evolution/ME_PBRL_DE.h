#include "Information_Bottleneck/IB_kernel.h"
#include "Discrete_Density_Evolution/Probability_Combination_Tool.h"
#include "Information_Bottleneck/stats.h"
#include <string>
#include <string.h>

class ME_PBRL_DE
{

// properties
public:
    std::vector<std::vector<double>> variable_node_degree;
    std::vector<std::vector<double>> check_node_degree;
    std::vector<std::vector<double>> vari_edge_deg_1, vari_edge_deg_2, vari_edge_deg_3;
    std::vector<std::vector<double>> check_edge_deg_1, check_edge_deg_2, check_edge_deg_3;
    std::vector<std::vector<double>> vari_representation_1, vari_representation_2;
    std::vector<std::vector<double>> check_representation_1, check_representation_2;
    std::vector<std::vector<double>> vari_threshold_1, vari_threshold_2;
    std::vector<std::vector<double>> check_threshold_1, check_threshold_2;
    std::vector<std::vector<double>> channel_threshold;
    std::vector<std::vector<double>> channel_representation;
    std::vector<std::vector<double>> vari_pmf_1, vari_pmf_2, vari_pmf_3;
    std::vector<std::vector<double>> check_pmf_1, check_pmf_2, check_pmf_3;
    std::vector<std::vector<double>> p_channel_pmf, np_channel_pmf;
    double sigma2;
    unsigned int max_iter;
    unsigned int quantization_size;
    double stop_threshold;
    double llr_combination_interval;
    unsigned ib_runtime;
    std::string pbrl_met_description;
    std::string suffix;

// methods
public:
    ME_PBRL_DE(std::string PBRL_MET_Description, unsigned Max_iter, unsigned Quansize, double Stop_treshold, double LLR_interval, unsigned IB_runtime, std::string Suffix);
    std::vector<std::vector<double>> calculate_output_distribution(std::vector<double> &distribution, const char type[]);
    void type_distribution_update(std::vector<std::vector<double>> edge_distribution, const char type[], int iter, int socket);
    bool read_decription();
};