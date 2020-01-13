#include "Discrete_Density_Evolution/ME_PBRL_DE.h"

ME_PBRL_DE::ME_PBRL_DE(const std::string PBRL_MET_Description, unsigned Max_iter,
                       unsigned Quansize, double Sigma2, double Stop_treshold,
                       double LLR_interval, unsigned IB_runtime, std::string Suffix)
{
    pbrl_met_description = PBRL_MET_Description;
    max_iter = Max_iter;
    quantization_size = Quansize;
    stop_threshold = Stop_treshold;
    llr_combination_interval = LLR_interval;
    ib_runtime = IB_runtime;
    suffix = Suffix;
    sigma2 = Sigma2;
}

std::vector<std::vector<double>> ME_PBRL_DE::calculate_output_distribution(std::vector<double> &distribution, const char type[])
{
    std::vector<std::vector<double>> output_joint(2);
    std::vector<std::vector<double>> first_input, second_input, combined_output, prob_llr_combined;
    std::vector<double> temp;
    if (strcmp(type, "vari") == 0)
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
           
            if (distribution[2] > 0)
            {
                second_input = check_pmf_1;
                for (int ii = 0; ii < int(distribution[2]); ii++)
                {
                    combined_output = prob_combination(first_input, second_input, type);
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input = prob_llr_combined;
                }
            }

            if (distribution[3] > 0)
            {
                second_input = check_pmf_2;
                for (int ii = 0; ii < int(distribution[3]); ii++)
                {
                    combined_output = prob_combination(first_input, second_input, type);
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input = prob_llr_combined;
                }
            }
            if (distribution[4] > 0)
            {
                second_input = check_pmf_3;
                for (int ii = 0; ii < int(distribution[4]); ii++)
                {
                    combined_output = prob_combination(first_input, second_input, type);
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input = prob_llr_combined;
                }
            }

            temp = distribution[0] * first_input[0];
            std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[0]));
            temp = distribution[0] * first_input[1];
            std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[1]));
            prob_sort(output_joint);
            std::cout<<"has been here 2 "<<std::endl;
            return output_joint;
        }
        else
        {
            std::cout << "Wrong Info: It's not a variable node distribution vector (length: " << distribution.size() << "), plz check again" << std::endl;
            return output_joint;
        }
    }
    else if (strcmp(type, "check") == 0)
    {
        std::cout << std::endl;
        if (distribution.size() == 4)
        {
            std::vector<unsigned> ind;
            std::vector<std::vector<double>> first_input_c1, first_input_c2, first_input_c3;
            //---------------------------caculate respectively---------------------------------
            //tyep 1
            if (distribution[1] != 0)
            {
                ind.push_back(0);
                first_input_c1 = vari_pmf_1;
                second_input = vari_pmf_1;
                for (int ii = 0; ii < int(distribution[1]) - 1; ii++)
                {
                    combined_output = prob_combination(first_input_c1, second_input, type);
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
                    combined_output = prob_combination(first_input_c2, second_input, type);
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
                    combined_output = prob_combination(first_input_c3, second_input, type);
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
                    combined_output = prob_combination(first_input_c1, first_input_c2, type);
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    return prob_llr_combined;
                }
                else if (ind[0] == 0 && ind[1] == 2)
                {
                    combined_output = prob_combination(first_input_c1, first_input_c3, type);
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    return prob_llr_combined;
                }
                else if (ind[0] == 1 && ind[1] == 2)
                {
                    combined_output = prob_combination(first_input_c2, first_input_c3, type);
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
        if (type =="vari")
            std::cout<<edge_distribution.size()<<"  "<<"has been here \n";
    }
    
    prob_sort(output_joint);
    llr_combined_output_joint = llr_combination(output_joint, llr_combination_interval);
    //std::cout << "Info: Type--"<<socket<<"--Node: "<<type<<"---before combined mi: " << it_mi(combined_vari_dist) << ";  after combined mi: " << it_mi(llr_combined_vari_dist) << std::endl;
    clipped_cvd = clip_prob(llr_combined_output_joint, pow(10, -10.0));
    ave_joinprob_llr(clipped_cvd, pow(10.0, -80.0));
    llr = llr_cal(clipped_cvd);
    llr_file_name = "vari_type_s" + std::to_string(socket) + "_llr_iteration_" + std::to_string(iter) + ".txt";
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

    if (strcmp(type, "vari") == 0)
    {
        switch (socket)
        {
        case 0:
            vari_representation_1.push_back(llr_cal(IB_ins.prob_join_xt));
            vari_threshold_1.push_back(IB_ins.threshold);
            vari_pmf_1 = IB_ins.prob_join_xt;
            std::cout << "Info: Iteration " << iter << ", socket " << socket << ", variable mutual information update: " << IB_ins.mi << std::endl;
            break;
        case 1:
            vari_representation_2.push_back(llr_cal(IB_ins.prob_join_xt));
            vari_threshold_2.push_back(IB_ins.threshold);
            vari_pmf_2 = IB_ins.prob_join_xt;
            std::cout << "Info: Iteration " << iter << ", socket " << socket << ", variable mutual information update: " << IB_ins.mi << std::endl;
            break;
        default:
            std::cout << "Wrong Info: In function (type_distribution_update), type is vari, but no socket " << socket << ".. please check again." << std::endl;
            break;
        }
    }
    else if (strcmp(type, "check") == 0)
    {
        switch (socket)
        {
        case 0:
            check_representation_1.push_back(llr_cal(IB_ins.prob_join_xt));
            check_threshold_1.push_back(IB_ins.threshold);
            check_pmf_1 = IB_ins.prob_join_xt;
            break;
        case 1:
            check_representation_2.push_back(llr_cal(IB_ins.prob_join_xt));
            check_threshold_2.push_back(IB_ins.threshold);
            check_pmf_2 = IB_ins.prob_join_xt;
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

int ME_PBRL_DE::density_evolution()
{
    //-----------------initial channel quantization--------------------------
    double most_left = -5.0;
    double most_right = 5.0;
    int partition_number = 3000;
    double small_offset = pow(10, -3);
    std::vector<std::vector<double>> temp_two_dim_vec;
    std::vector<std::vector<double>> channel_observation = gaussian_disretization(most_left, most_right, partition_number, sigma2);
    IB_kernel channel_IB(channel_observation, quantization_size, ib_runtime);
    channel_IB.Progressive_MMI();
    np_channel_pmf = channel_IB.prob_join_xt;
    p_channel_pmf = {{0.5 * 1 / (1 + exp(small_offset)), 0.5 * exp(small_offset) / (1 + exp(small_offset))},
                     {0.5 * exp(small_offset) / (1 + exp(small_offset)), 0.5 * 1 / (1 + exp(small_offset))}};

    //-----------------initial variable node edge distribution---------------
    std::vector<double> channel_dist(2, 0); // (>0)->(not puncture)->(0) || (<0)->punc->(1)
    for (unsigned ii = 0; ii < vari_edge_deg_1.size(); ii++)
    {
        if (vari_edge_deg_1[ii][1] > 0)
            channel_dist[0] += vari_edge_deg_1[ii][0];
        else
            channel_dist[1] += vari_edge_deg_1[ii][0];
    }
    if (channel_dist[0] + channel_dist[1] < 1 - pow(10, -7))
    {
        std::cout << "Wrong Info: Summation of variable edge dstirbution of socket 1 is not equal to 1 rather " << channel_dist[0] + channel_dist[1] << ", plz check and try again..." << std::endl;
        return -1;
    }
    temp_two_dim_vec.clear();
    temp_two_dim_vec.push_back(channel_dist[0] * channel_observation[0]);
    temp_two_dim_vec.push_back(channel_dist[0] * channel_observation[1]);
    temp_two_dim_vec[0].push_back(channel_dist[1] * 0.5 * 1 / (1 + exp(small_offset)));
    temp_two_dim_vec[0].push_back(channel_dist[1] * 0.5 * exp(small_offset) / (1 + exp(small_offset)));
    temp_two_dim_vec[1].push_back(channel_dist[1] * 0.5 * exp(small_offset) / (1 + exp(small_offset)));
    temp_two_dim_vec[1].push_back(channel_dist[1] * 0.5 * 1 / (1 + exp(small_offset)));
    prob_sort(temp_two_dim_vec);
    IB_kernel channel_IB_1(temp_two_dim_vec, quantization_size, ib_runtime);
    channel_IB_1.Progressive_MMI();
    vari_pmf_1 = channel_IB_1.prob_join_xt;
    channel_dist.assign(2, 0);
    for (unsigned ii = 0; ii < vari_edge_deg_2.size(); ii++)
    {
        if (vari_edge_deg_2[ii][1] > 0)
            channel_dist[0] += vari_edge_deg_2[ii][0];
        else
            channel_dist[1] += vari_edge_deg_2[ii][0];
    }
    if (channel_dist[0] + channel_dist[1] < 1 - pow(10, -7))
    {
        std::cout << "Wrong Info: Summation of variable edge dstirbution of socket 2 is not equal to 1 rather " << channel_dist[0] + channel_dist[1] << ", plz check and try again..." << std::endl;
        return -1;
    }
    temp_two_dim_vec.clear();
    temp_two_dim_vec.push_back(channel_dist[0] * channel_observation[0]);
    temp_two_dim_vec.push_back(channel_dist[0] * channel_observation[1]);
    temp_two_dim_vec[0].push_back(channel_dist[1] * 0.5 * 1 / (1 + exp(small_offset)));
    temp_two_dim_vec[0].push_back(channel_dist[1] * 0.5 * exp(small_offset) / (1 + exp(small_offset)));
    temp_two_dim_vec[1].push_back(channel_dist[1] * 0.5 * exp(small_offset) / (1 + exp(small_offset)));
    temp_two_dim_vec[1].push_back(channel_dist[1] * 0.5 * 1 / (1 + exp(small_offset)));
    prob_sort(temp_two_dim_vec);
    IB_kernel channel_IB_2(temp_two_dim_vec, quantization_size, ib_runtime);
    channel_IB_2.Progressive_MMI();
    vari_pmf_2 = channel_IB_2.prob_join_xt;
    vari_pmf_3 = np_channel_pmf;

    channel_threshold.push_back(channel_IB_1.threshold);
    channel_threshold.push_back(channel_IB_2.threshold);
    channel_threshold.push_back(channel_IB.threshold);
    channel_representation.push_back(llr_cal(channel_IB_1.prob_join_xt));
    channel_representation.push_back(llr_cal(channel_IB_2.prob_join_xt));
    channel_representation.push_back(llr_cal(channel_IB.prob_join_xt));
    std::cout << "Info: Finished channel quantization ..." << std::endl;
    //---------------------density evolution------------------------------------
    for (unsigned ii = 0; ii < max_iter; ii++)
    {
        //-----------------update check edge pmf---------------------------------
        type_distribution_update(check_edge_deg_1, "check", ii, 0);
        type_distribution_update(check_edge_deg_2, "check", ii, 1);
        //-----------------update variable node pmf------------------------------
        type_distribution_update(vari_edge_deg_1, "vari", ii, 0);
        type_distribution_update(vari_edge_deg_2, "vari", ii, 1);
    }

    RQF_output();
    return 1;
}

void ME_PBRL_DE::RQF_output()
{
    //-------------------First Writeout channel_information---------------------
    std::string channelquan_file = "channel_quantizer_" + suffix + ".txt";
    std::string channelrec_file = "channel_reconstruction_" + suffix + ".txt";
    std::ofstream handle_channel_quantizer(channelquan_file);
    std::vector<std::vector<double>> cur_two_dim_vec;
    if (handle_channel_quantizer.is_open())
    {
        handle_channel_quantizer << quantization_size << "  " << std::endl;
        cur_two_dim_vec = channel_threshold;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                handle_channel_quantizer << bb << "  ";
            handle_channel_quantizer << std::endl;
        }
        handle_channel_quantizer.close();
    }
    else
    {
        std::cout << "Fail to write channel quantizer ..." << std::endl;
    }

    std::ofstream handle_channel_recons(channelrec_file);
    if (handle_channel_recons.is_open())
    {
        handle_channel_recons << quantization_size << "  " << std::endl;
        cur_two_dim_vec = channel_representation;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                handle_channel_recons << bb << "  ";
            handle_channel_recons << std::endl;
        }
        handle_channel_recons.close();
    }
    else
    {
        std::cout << "Fail to write channel reconstruction file ..." << std::endl;
    }

    std::ofstream outpufile, handle_reconstrction;
    std::string threholdfile, recons_file;

    //------------------output socket 1 information---------------------------------
    threholdfile = "threshold_s1" + suffix + ".txt";
    recons_file = "reconstruction_s1" + suffix + ".txt";
    outpufile.open(threholdfile);
    if (outpufile.is_open())
    {
        outpufile << quantization_size << " " << max_iter << "  " << std::endl;
        cur_two_dim_vec = check_threshold_1;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                outpufile << bb << "  ";
            outpufile << std::endl;
        }
        cur_two_dim_vec = vari_threshold_1;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                outpufile << bb << "  ";
            outpufile << std::endl;
        }

        outpufile.close();
    }

    //Part 2: Reconstruction function
    handle_reconstrction.open(recons_file);
    if (handle_reconstrction.is_open())
    {
        handle_reconstrction << quantization_size << "  " << max_iter << "  " << std::endl;
        cur_two_dim_vec = check_representation_1;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                handle_reconstrction << bb << "  ";
            handle_reconstrction << std::endl;
        }
        cur_two_dim_vec = vari_representation_1;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                handle_reconstrction << bb << "  ";
            handle_reconstrction << std::endl;
        }
        handle_reconstrction.close();
    }
    else
    {
        std::cout << "Reconstruction failed to wirte ..." << std::endl;
    }

    //--------------------output sockect two information-----------------------------------
    threholdfile = "threshold_s2" + suffix + ".txt";
    recons_file = "reconstruction_s2" + suffix + ".txt";
    outpufile.open(threholdfile);
    if (outpufile.is_open())
    {
        outpufile << quantization_size << "  " << max_iter << "  " << std::endl;
        cur_two_dim_vec = check_threshold_2;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                outpufile << bb << "  ";
            outpufile << std::endl;
        }
        cur_two_dim_vec = vari_threshold_2;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                outpufile << bb << "  ";
            outpufile << std::endl;
        }

        outpufile.close();
    }

    //Part 2: Reconstruction function
    handle_reconstrction.open(recons_file);
    if (handle_reconstrction.is_open())
    {
        handle_reconstrction << quantization_size << "  " << max_iter << "  " << std::endl;
        cur_two_dim_vec = check_representation_2;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                handle_reconstrction << bb << "  ";
            handle_reconstrction << std::endl;
        }
        cur_two_dim_vec = vari_representation_2;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                handle_reconstrction << bb << "  ";
            handle_reconstrction << std::endl;
        }
        handle_reconstrction.close();
    }
    else
    {
        std::cout << "Reconstruction failed to wirte ..." << std::endl;
    }
}