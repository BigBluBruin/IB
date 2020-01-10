#include "Discrete_Density_Evolution/ME_PBRL_DE.h"






std::vector<std::vector<double>> ME_PBRL_DE::calculate_output_distribution(std::vector<double> & distribution, const char type[])
{
    std::vector<std::vector<double>> output_joint(2);
    std::vector<std::vector<double>> first_input,second_input,combined_output,prob_llr_combined;
    std::vector<double> temp;

    if (type == "vari")
    {
        if (distribution.size() == 5)
        {
            if (distribution[1] > 0) //>0 not pucture
            {
                first_input = np_channel_pmf;
            }
            else //<0 puncture
            {
                 first_input = {{0.25,0.25},{0.25,0.25}};
            }
            if (distribution[2] != 0)
            {
                for (unsigned ii = 0; ii < distribution[2]; ii++)
                {
                    second_input = check_pmf_1;
                    combined_output = prob_combination(first_input, second_input, "check");
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input = prob_llr_combined;
                }
                temp = distribution[0] * first_input[0];
                std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[0]));
                temp = distribution[0] * first_input[1];
                std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[1]));
            }
            if (distribution[3] != 0)
            {
                for (unsigned ii = 0; ii < distribution[3]; ii++)
                {
                    second_input = check_pmf_2;
                    combined_output = prob_combination(first_input, second_input, "check");
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input = prob_llr_combined;
                }
                temp = distribution[0] * first_input[0];
                std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[0]));
                temp = distribution[0] * first_input[1];
                std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[1]));
            }
            if (distribution[4] != 0)
            {
                for (unsigned ii = 0; ii < distribution[4]; ii++)
                {
                    second_input = check_pmf_2;
                    combined_output = prob_combination(first_input, second_input, "check");
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input = prob_llr_combined;
                }
                temp = distribution[0] * first_input[0];
                std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[0]));
                temp = distribution[0] * first_input[1];
                std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[1]));
            }
        }
        else
        {
            std::cout<<"Wrong Info: It's not a variable node distribution vector (length: "<<distribution.size()<<"), plz check again"<<std::endl;
        }
        

    }
    else if (type=="check")
    {
        if(distribution.size()==4)
        {

        }
        else
        {
            std::cout<<"Wrong Info: It's not a variable node distribution vector (length: "<<distribution.size()<<"), plz check again"<<std::endl;
        }
        
    }
    else
    {
        std::cout<<"Neiter vari nor check, plz check again..."<<std::endl;
    }
    
    
}