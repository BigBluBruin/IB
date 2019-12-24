#include "Quantize_Continuous_DE.h"

bool QCDE(const std::string filename)
{
    std::ifstream handle_file(filename);
    if (handle_file.is_open())
    {
        unsigned total_points, max_iter;
        handle_file >> total_points >> max_iter;
        std::vector<double> llr_val(total_points);
        for (unsigned index1 = 0; index1 < total_points; index1++)
        {
            handle_file >> llr_val[index1];
        }
        std::vector<std::vector<double>> check_llr(max_iter);
        for (unsigned index1 = 0; index1 < max_iter; index1++)
        {
            check_llr[index1].assign(total_points, -1);
            for (unsigned index2 = 0; index2 < total_points; index2++)
            {
                handle_file >> check_llr[index1][index2];
            }
        }
        std::vector<std::vector<double>> vari_llr(max_iter);
        for (unsigned index1 = 0; index1 < max_iter; index1++)
        {
            vari_llr[index1].assign(total_points, -1);
            for (unsigned index2 = 0; index2 < total_points; index2++)
            {
                handle_file >> vari_llr[index1][index2];
            }
        }
        std::vector<std::vector<double>> joint_prob_xy, joint_prob_xt;
        std::vector<std::vector<double>> check_threshold,vari_threshold;
        std::vector<std::vector<double>> check_recons,vari_recons;
    }
    else
    {
    }
}