#include "Discrete_Density_Evolution/ME_PBRL_DE.h"

ME_PBRL_DE::ME_PBRL_DE(std::string PBRL_MET_Description = " ",
                       unsigned Max_iter = 50, unsigned Quansize = 16, double Stop_treshold = 0, double LLR_interval = 0.01, unsigned IB_runtime = 50, std::string Suffix = " ")
{
    pbrl_met_description=PBRL_MET_Description;
    max_iter = Max_iter;
    quantization_size = Quansize;
    stop_threshold = Stop_treshold;
    llr_combination_interval = LLR_interval;
    ib_runtime = IB_runtime;
    suffix = Suffix;
}

std::vector<std::vector<double>> ME_PBRL_DE::calculate_output_distribution(std::vector<double> &distribution, const char type[])
{
    std::vector<std::vector<double>> output_joint(2);
    std::vector<std::vector<double>> first_input, second_input, combined_output, prob_llr_combined;
    std::vector<double> temp;

    if (strcmp(type,"vari")==0)
    {
        if (distribution.size() == 5)
        {
            if (distribution[1] > 0) //>0 not pucture
            {
                first_input = np_channel_pmf;
            }
            else //<0 puncture
            {
                first_input = {{0.5 / (1 + exp(0.001)), 0.5 * exp(0.001) / (1 + exp(0.001))}, {0.5 * exp(0.001) / (1 + exp(0.001)), 0.5 / (1 + exp(0.001))}}; //a 0.001dB offset
            }
            if (distribution[2] != 0)
            {
                second_input = check_pmf_1;
                for (int ii = 0; ii < int(distribution[2]); ii++)
                {
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
                second_input = check_pmf_2;
                for (int ii = 0; ii < int(distribution[3]); ii++)
                {
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
                second_input = check_pmf_2;
                for (int ii = 0; ii < int(distribution[4]); ii++)
                {
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
            prob_sort(output_joint);
            return output_joint;
        }
        else
        {
            std::cout << "Wrong Info: It's not a variable node distribution vector (length: " << distribution.size() << "), plz check again" << std::endl;
            return output_joint;
        }
    }
    else if (strcmp(type,"check")==0)
    {
        if (distribution.size() == 4)
        {
            std::vector<unsigned> ind;
            std::vector<std::vector<double>> first_input_c1(2), first_input_c2(2), first_input_c3(2);
            //---------------------------caculate respectively---------------------------------
            //tyep 1
            if (distribution[1] != 0)
            {
                ind.push_back(0);
                first_input_c1 = vari_pmf_1;
                second_input = vari_pmf_1;
                for (int ii = 0; ii < int(distribution[1]) - 1; ii++)
                {
                    combined_output = prob_combination(first_input_c1, second_input, "check");
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input_c1 = prob_llr_combined;
                }
            }
            //type 2
            if (distribution[2] != 0)
            {
                ind.push_back(1);
                first_input_c2 = vari_pmf_2;
                second_input = vari_pmf_2;
                for (int ii = 0; ii < int(distribution[2]) - 1; ii++)
                {
                    combined_output = prob_combination(first_input_c2, second_input, "check");
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input_c2 = prob_llr_combined;
                }
            }
            //type 3
            if (distribution[3] != 0)
            {
                ind.push_back(2);
                first_input_c3 = vari_pmf_3;
                second_input = vari_pmf_3;
                for (int ii = 0; ii < int(distribution[3]) - 1; ii++)
                {
                    combined_output = prob_combination(first_input_c3, second_input, "check");
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input_c3 = prob_llr_combined;
                }
            }
            //-----------------------------------------------------------------------------------------

            //-----------------------------fianl calculation-------------------------------------------
            if (ind.size() == 1)
            {
                if (ind[0] == 0)
                {
                    return first_input_c1;
                }
                else if (ind[0] == 1)
                {
                    return first_input_c2;
                }
                else if (ind[0] == 2)
                {
                    return first_input_c3;
                }
                else
                {
                    std::cout << "Wrong Info: No such one check node, ind is: " << ind[0] << ". plz check..." << std::endl;
                    return output_joint;
                }
            }
            else if (ind.size() == 2)
            {
                if (ind[0] == 0 && ind[1] == 1)
                {
                    combined_output = prob_combination(first_input_c1, first_input_c2, "check");
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    return prob_llr_combined;
                }
                else if (ind[0] == 0 && ind[1] == 2)
                {
                    combined_output = prob_combination(first_input_c1, first_input_c3, "check");
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    return prob_llr_combined;
                }
                else if (ind[0] == 1 && ind[1] == 2)
                {
                    combined_output = prob_combination(first_input_c2, first_input_c3, "check");
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    return prob_llr_combined;
                }
                else
                {
                    std::cout << "Wrong Info: No such two check node, inds are: (" << ind[0] << ", " << ind[1] << "). plz check..." << std::endl;
                    return output_joint;
                }
            }
            else if (ind.size() == 3)
            {
                std::cout << "Wrong Info: No 3 combinations can happen in this case, plz check..." << std::endl;
                return output_joint;
            }
            else
            {
                std::cout << "Wrong Info: No such check node combinations, ind size is: " << ind.size() << ". plz check..." << std::endl;
                return output_joint;
            }

            //-----------------------------------------------------------------------------------------
        }
        else
        {
            std::cout << "Wrong Info: It's not a check node distribution vector (length: " << distribution.size() << ", no equal to 4), plz check again" << std::endl;
            return output_joint;
        }
    }
    else
    {
        std::cout << "Neiter vari nor check, plz check again..." << std::endl;
        return output_joint;
    }
}

void ME_PBRL_DE::type_distribution_update(std::vector<std::vector<double>> edge_distribution,
                                                                      const char type[], int iter, int socket)
{
    std::vector<std::vector<double>> output_joint(2), cur_output;
    std::vector<std::vector<double>> llr_combined_output_joint, clipped_cvd;
    std::vector<double> llr;
    std::ofstream llr_file;
    std::string llr_file_name;
    for (unsigned ii = 0; ii < edge_distribution.size(); ii++)
    {
        cur_output = calculate_output_distribution(edge_distribution[ii], type);
        std::copy(cur_output[0].begin(), cur_output[0].end(), std::back_inserter(output_joint[0]));
        std::copy(cur_output[1].begin(), cur_output[1].end(), std::back_inserter(output_joint[1]));
    }
    prob_sort(output_joint);
    llr_combined_output_joint = llr_combination(output_joint, llr_combination_interval);
    //std::cout << "Info: Type--"<<socket<<"--Node: "<<type<<"---before combined mi: " << it_mi(combined_vari_dist) << ";  after combined mi: " << it_mi(llr_combined_vari_dist) << std::endl;
    clipped_cvd = clip_prob(llr_combined_output_joint, pow(10, -10.0));
    ave_joinprob_llr(clipped_cvd, pow(10.0, -80.0));
    llr = llr_cal(clipped_cvd);
    llr_file_name = "vari_type_" + std::to_string(socket) + "_llr_iteration_" + std::to_string(iter) + ".txt";
    llr_file.open(llr_file_name);
    if (llr_file.is_open())
    {
        for (unsigned index = 0; index < llr.size(); index++)
        {
            llr_file << llr[index] << "  " << clipped_cvd[0][index] + clipped_cvd[1][index] << std::endl;
        }
    }
    llr_file.close();
    IB_kernel IB_ins(clipped_cvd, quantization_size, ib_runtime);
    IB_ins.Progressive_MMI();

    if (strcmp(type,"vari")==0)
    {
        switch (socket)
        {
        case 0:
            vari_representation_1.push_back(llr_cal(IB_ins.prob_join_xt));
            vari_threshold_1.push_back(IB_ins.threshold);
            vari_pmf_1=IB_ins.prob_join_xt;
            break;
        case 1:
            vari_representation_2.push_back(llr_cal(IB_ins.prob_join_xt));
            vari_threshold_2.push_back(IB_ins.threshold);
            vari_pmf_2=IB_ins.prob_join_xt;
            break;
        default:
            std::cout << "Wrong Info: In function (type_distribution_update), type is vari, but no socket " << socket << ".. please check again." << std::endl;
            break;
        }
    }
    else if (strcmp(type,"check")==0)
    {
        switch (socket)
        {
        case 0:
            check_representation_1.push_back(llr_cal(IB_ins.prob_join_xt));
            check_threshold_1.push_back(IB_ins.threshold);
            check_pmf_1=IB_ins.prob_join_xt;
            break;
        case 1:
            check_representation_2.push_back(llr_cal(IB_ins.prob_join_xt));
            check_threshold_2.push_back(IB_ins.threshold);
            vari_pmf_2=IB_ins.prob_join_xt;
            break;
        default:
            std::cout << "Wrong Info: In function (type_distribution_update), type is check, but no socket " << socket << ".. please check again." << std::endl;
            break;
        }
    }
    else
    {
        std::cout << "Wrong Info: In function (type_distribution_update), type is neither vari or check, it is " << type << ". Plz check again." << std::endl;
    }
}

bool ME_PBRL_DE::read_decription()
{
    std::ifstream myfile(pbrl_met_description);
    if (myfile.is_open())
    {
        int cardi; //number of terms in each distribution

        //-------variable node -------------------------
        myfile >> cardi;
        variable_node_degree.resize(cardi, std::vector<double>(5, -1));
        for (unsigned ii = 0; ii < variable_node_degree.size(); ii++)
        {
            for (unsigned jj = 0; jj < 5; jj++)
            {
                myfile >> variable_node_degree[ii][jj];
            }
        }

        //------check node----------------------------
        myfile >> cardi;
        check_node_degree.resize(cardi, std::vector<double>(4, -1));
        for (unsigned ii = 0; ii < check_node_degree.size(); ii++)
        {
            for (unsigned jj = 0; jj < 4; jj++)
            {
                myfile >> check_node_degree[ii][jj];
            }
        }

        //------variable edge socket 1----------------------
        myfile >> cardi;
        vari_edge_deg_1.resize(cardi, std::vector<double>(5, -1));
        for (unsigned ii = 0; ii < vari_edge_deg_1.size(); ii++)
        {
            for (unsigned jj = 0; jj < 5; jj++)
            {
                myfile >> vari_edge_deg_1[ii][jj];
            }
        }

        //------variable edge  socket 2-----------------------
        myfile >> cardi;
        vari_edge_deg_2.resize(cardi, std::vector<double>(5, -1));
        for (unsigned ii = 0; ii < vari_edge_deg_2.size(); ii++)
        {
            for (unsigned jj = 0; jj < 5; jj++)
            {
                myfile >> vari_edge_deg_2[ii][jj];
            }
        }

        //------variable edge  socket 3-----------------------
        myfile >> cardi;
        vari_edge_deg_3.resize(cardi, std::vector<double>(5, -1));
        for (unsigned ii = 0; ii < vari_edge_deg_3.size(); ii++)
        {
            for (unsigned jj = 0; jj < 5; jj++)
            {
                myfile >> vari_edge_deg_3[ii][jj];
            }
        }

        //-------check edge socket 1---------------------------
        myfile >> cardi;
        check_edge_deg_1.resize(cardi, std::vector<double>(4, -1));
        for (unsigned ii = 0; ii < check_edge_deg_1.size(); ii++)
        {
            for (unsigned jj = 0; jj < 4; jj++)
            {
                myfile >> check_edge_deg_1[ii][jj];
            }
        }

        //-------check edge socket 2---------------------------
        myfile >> cardi;
        check_edge_deg_2.resize(cardi, std::vector<double>(4, -1));
        for (unsigned ii = 0; ii < check_edge_deg_2.size(); ii++)
        {
            for (unsigned jj = 0; jj < 4; jj++)
            {
                myfile >> check_edge_deg_2[ii][jj];
            }
        }

        //-------check edge socket 3---------------------------
        myfile >> cardi;
        check_edge_deg_3.resize(cardi, std::vector<double>(4, -1));
        for (unsigned ii = 0; ii < check_edge_deg_3.size(); ii++)
        {
            for (unsigned jj = 0; jj < 4; jj++)
            {
                myfile >> check_edge_deg_3[ii][jj];
            }
        }

        myfile.close();
        std::cout << "Info: Load information from " << pbrl_met_description << " successfully..." << std::endl;
        return true;
    }
    else
    {
        std::cout << "Wrong Info: Can't find the file: " << pbrl_met_description << ". Plz check file again..." << std::endl;
        return false;
    }
}
