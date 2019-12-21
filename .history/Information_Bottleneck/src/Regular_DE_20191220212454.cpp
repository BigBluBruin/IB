#include "Discrete_Density_Evolution/Regular_DE.h"

Regular_DE::Regular_DE(unsigned int Dc=0, unsigned int Dv=0, double Sigma2=0.0, unsigned int Max_iter=0, unsigned int Quantization_size=0, double Threshold=0)
{
    dc=Dc;
    dv=Dv;
    sigma2=Sigma2;
    max_iter=Max_iter;
    quantization_size=Quantization_size;
    threshold=Threshold;

    
}

int Regular_DE::Discrete_Density_Evolution()
{
    std::vector<std::vector<double>> channel_observation = gaussian_disretization(-5, 5, 3000, sigma2);
    IB_kernel channel_IB(channel_observation, quantization_size,500);
    channel_IB.smIB();
    std::vector<std::vector<double>> first_input = channel_IB.prob_join_xt;
    std::vector<std::vector<double>> second_input = first_input;
    std::vector<std::vector<double>> combined_output, prob_llr_combined;
    std::cout<<"finished channel quantization...."<<std::endl;
    for (unsigned iter = 0; iter < max_iter; iter++)
    {
        //---------------check node--------------------
        for (unsigned ii = 0; ii < dc - 2; ii++)
        {
            combined_output = prob_combination(first_input, second_input, "check");
            //---------------check part-----------------------------
            //------------------------------------------------------
            prob_sort(combined_output);
            prob_llr_combined = llr_combination(combined_output, pow(10,-4));
            //ave_joinprob(prob_llr_combined);
            ave_joinprob_llr( prob_llr_combined, pow(10.0,-80.0));
            first_input = prob_llr_combined;
        }
        IB_kernel check_IB(prob_llr_combined, quantization_size,2000);
        check_IB.smIB();
        check_representation.push_back(llr_cal(check_IB.prob_join_xt));
        check_threshold.push_back(check_IB.threshold);
        //--------------variable node-------------------
        first_input = channel_IB.prob_join_xt;
        second_input = check_IB.prob_join_xt;
        for (unsigned ii = 0; ii < dv - 1; ii++)
        {
            combined_output = prob_combination(first_input, second_input, "vari");
            prob_sort(combined_output);
            prob_llr_combined = llr_combination(combined_output, 0);
            //ave_joinprob(prob_llr_combined);
            ave_joinprob_llr( prob_llr_combined, pow(10.0,-80.0));
            first_input = prob_llr_combined;
        }
        std::vector<double> llr=llr_cal(first_input);
        std::string llr_file_name="llr_iteration"+std::to_string(iter) +".txt";
        std::ofstream llr_file(llr_file_name);
        if(llr_file.is_open())
        {
            for(unsigned index=0; index<llr.size();index++)
            {
                llr_file<<llr[index]<<"  "<<first_input[0][index]+first_input[1][index]<<std::endl;
            }
        }
        llr_file.close();
        IB_kernel vari_IB(prob_llr_combined, quantization_size,2000);
        vari_IB.smIB();
        vari_representation.push_back(llr_cal(vari_IB.prob_join_xt));
        vari_threshold.push_back(vari_IB.threshold);
        //---------------check node------------------------
        first_input=vari_IB.prob_join_xt;
        second_input=first_input;
        //--------------Output Info------------------------
        std::cout<<"DE Info:  iter--"<<iter+1<<"--mi--"<<vari_IB.mi<<std::endl;
        std::cout<<"Variable node threshold:"<<std::endl;
        for(const auto & aa: vari_IB.threshold)
            std::cout << aa << "  ";
        std::cout << std::endl;
        if (1 - vari_IB.mi < threshold)
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
