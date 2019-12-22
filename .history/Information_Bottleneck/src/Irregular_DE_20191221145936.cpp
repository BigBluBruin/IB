#include "Discrete_Density_Evolution/Irregular_DE.h"

Irregular_DE::Irregular_DE(std::vector<double> Check_edge_dist={0.0}, std::vector<double> Vari_edge_dist={0.0}, double Sigma2=0.0, unsigned int Max_iter=0,unsigned int Quantization_size=0,double Stop_threshold=0.0,double Llr_combination_interval=0.0,unsigned Ib_runtime=0.0)
{
    check_edge_dist=Check_edge_dist;
    vari_edge_dist=Vari_edge_dist;
    sigma2=Sigma2;
    max_iter=Max_iter;
    quantization_size=Quantization_size;
    stop_threshold=Stop_threshold;
    llr_combination_interval=Llr_combination_interval;
    ib_runtime=Ib_runtime;
}


int Irregular_DE::Discrete_Density_Evolution()
{
    std::vector<std::vector<double>> combined_check_dist, combined_vari_dist,llr_combined_check_dist, llr_combined_vari_dist;
    std::vector<double> temp;
    std::vector<std::vector<double>> channel_observation = gaussian_disretization(-5, 5, 3000, sigma2);
    IB_kernel channel_IB(channel_observation, quantization_size, 2000);
    channel_IB.smIB();
    std::vector<std::vector<double>> first_input = channel_IB.prob_join_xt;
    std::vector<std::vector<double>> second_input = first_input;
    std::vector<std::vector<double>> combined_output, prob_llr_combined;
    std::cout<<"finished channel quantization...."<<std::endl;

    //-------------------start iteration------------------------------------------
    for (unsigned iter = 0; iter < max_iter; iter++)
    {
        combined_check_dist.clear();
        combined_vari_dist.clear();


        //---------------check node--------------------
        for (unsigned ii = 0; ii < check_edge_dist.size() - 2; ii++)
        {
            combined_output = prob_combination(first_input, second_input, "check");
            prob_sort(combined_output);
            prob_llr_combined = llr_combination(combined_output, 0.001);
            ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
            first_input = prob_llr_combined;
            if(check_edge_dist[ii+2]!=0)
            {
                temp=check_edge_dist[ii+2]*first_input[0];
                std::copy(temp.begin(),temp.end(),std::back_inserter(combined_check_dist[0]));
                temp=check_edge_dist[ii+2]*first_input[1];
                std::copy(temp.begin(),temp.end(),std::back_inserter(combined_check_dist[1]));
            }
        }
        prob_sort(combined_check_dist);
        llr_combined_check_dist = llr_combination(combined_check_dist, 0.001);
        std::vector<double> llr = llr_cal(llr_combined_check_dist);
        std::string llr_file_name = "check_llr_iteration_" + std::to_string(iter) + ".txt";
        std::ofstream llr_file(llr_file_name);
        if (llr_file.is_open())
        {
            for (unsigned index = 0; index < llr.size(); index++)
            {
                llr_file << llr[index] << "  " << llr_combined_check_dist[0][index] + llr_combined_check_dist[1][index] << std::endl;
            }
        }
        llr_file.close();
        IB_kernel check_IB(llr_combined_check_dist, quantization_size, 2000);
        check_IB.smIB();
        check_representation.push_back(llr_cal(check_IB.prob_join_xt));
        check_threshold.push_back(check_IB.threshold);
        //--------------------------------------------------------------


        //--------------variable node-------------------
        first_input = channel_IB.prob_join_xt;
        second_input = check_IB.prob_join_xt;
        for (unsigned ii = 0; ii < vari_edge_dist.size() - 1; ii++)
        {
            combined_output = prob_combination(first_input, second_input, "vari");
            prob_sort(combined_output);
            prob_llr_combined = llr_combination(combined_output, 0.001);
            ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
            first_input = prob_llr_combined;
            if (vari_edge_dist[ii + 1] != 0)
            {
                temp = vari_edge_dist[ii + 1] * first_input[0];
                std::copy(temp.begin(), temp.end(), std::back_inserter(combined_vari_dist[0]));
                temp = vari_edge_dist[ii + 1] * first_input[1];
                std::copy(temp.begin(), temp.end(), std::back_inserter(combined_vari_dist[1]));
            }
        }
        prob_sort(combined_vari_dist);
        llr_combined_vari_dist = llr_combination(combined_vari_dist, 0.001);
        llr = llr_cal(llr_combined_vari_dist);
        llr_file_name = "vari_llr_iteration_" + std::to_string(iter) + ".txt";
        llr_file.open(llr_file_name);
        if (llr_file.is_open())
        {
            for (unsigned index = 0; index < llr.size(); index++)
            {
                llr_file << llr[index] << "  " << llr_combined_vari_dist[0][index] + llr_combined_vari_dist[1][index] << std::endl;
            }
        }
        llr_file.close();
        IB_kernel vari_IB(llr_combined_vari_dist, quantization_size, 2000);
        vari_IB.smIB();
        vari_representation.push_back(llr_cal(vari_IB.prob_join_xt));
        vari_threshold.push_back(vari_IB.threshold);
        //------------------------------------------------------------------------

        //---------------check node------------------------
        first_input = vari_IB.prob_join_xt;
        second_input = first_input;
        //--------------Output Info------------------------
        std::cout << "DE Info:  iter--" << iter + 1 << "--mi--" << vari_IB.mi << std::endl;
        std::cout << "Variable node threshold:" << std::endl;
        for (const auto &aa : vari_IB.threshold)
            std::cout << aa << "  ";
        std::cout << std::endl;
        if (1 - vari_IB.mi < stop_threshold)
        {
            if (iter < max_iter - 3)
            {
                std::cout << "converge to fast..." << std::endl;
                return 1;
            }
            else
            {
                if (iter == max_iter - 1)
                {
                    std::cout << "find threshold" << std::endl;
                    std::ofstream outpufile("threshold.txt");
                    if (outpufile.is_open())
                    {
                        for (const auto &iter : check_threshold)
                        {
                            for (const auto &term : iter)
                                outpufile << term << "  ";
                            outpufile << std::endl;
                        }
                        outpufile << std::endl
                                  << std::endl;
                        for (const auto &iter : vari_threshold)
                        {
                            for (const auto &term : iter)
                                outpufile << term << "  ";
                            outpufile << std::endl;
                        }
                    }
                    return 3;
                }
            }
        }
        else
        {
            if (iter == max_iter - 1)
            {

                std::cout << "converge to slow..." << std::endl;
                return 2;
            }
        }
    }
}