#include "Discrete_Density_Evolution/ME2_PBRL_DE.h"

ME2_PBRL_DE::ME2_PBRL_DE(const std::string PBRL_MET_Description, unsigned Max_iter,
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

std::vector<std::vector<double>> ME2_PBRL_DE::calculate_output_distribution(std::vector<double> &distribution, const char type[])
{
    std::vector<std::vector<double>> output_joint(2);
    std::vector<std::vector<double>> first_input, second_input, combined_output, prob_llr_combined;
    std::vector<double> temp;
    if (strcmp(type, "vari") == 0)
    {
        if (distribution.size() == 4)
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
            temp = distribution[0] * first_input[0];
            std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[0]));
            temp = distribution[0] * first_input[1];
            std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[1]));
            prob_sort(output_joint);
            return output_joint;
        }
        else
        {
            std::cout << "Wrong Info: It's not a variable node distribution vector (length: " << distribution.size() << "), plz check again" << std::endl;
            return output_joint;
        }
    }
    else if ((strcmp(type, "check") == 0))
    {
        //std::cout<<distribution[0]<<std::endl;
        if (distribution.size() == 3)
        {
            std::vector<unsigned> ind;
            std::vector<std::vector<double>> first_input_c1, first_input_c2;
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
                    std::cout<<"Wrong Info: Check node distribution from second socked as at most 1, here the number is "<< int(distribution[2])<<". Please check again..."<<std::endl;
                    combined_output = prob_combination(first_input_c2, second_input, type);
                    prob_sort(combined_output);
                    prob_llr_combined = llr_combination(combined_output, 0.001);
                    ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                    first_input_c2 = prob_llr_combined;
                }
            }
            //-----------------------------------------------------------------------------------------
            
            //-----------------------------fianl calculation-------------------------------------------
            if (ind.size() == 1)
            {
                if (ind[0] == 0)
                {
                   
                    return distribution[0]*first_input_c1;
                }
                else if (ind[0] == 1)
                {
                    return distribution[0]*first_input_c2;
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
                    return distribution[0]*prob_llr_combined;
                }
                else
                {
                    std::cout << "Wrong Info: No such two check node, inds are: (" << ind[0] << ", " << ind[1] << "). plz check..." << std::endl;
                    return output_joint;
                }
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
    else if ((strcmp(type, "check_min") == 0))
    {
        //std::cout<<distribution[0]<<std::endl;
        if (distribution.size() == 3)
        {
            std::vector<unsigned> ind;
            std::vector<std::vector<double>> first_input_c1, first_input_c2;
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
                    first_input_c1 = combined_output;
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
                    std::cout << "Wrong Info: Check node distribution from second socked as at most 1, here the number is " << int(distribution[2]) << ". Please check again..." << std::endl;
                    combined_output = prob_combination(first_input_c2, second_input, type);
                    prob_sort(combined_output);
                    first_input_c2 = combined_output;
                }
            }
            //-----------------------------------------------------------------------------------------

            //-----------------------------fianl calculation-------------------------------------------
            if (ind.size() == 1)
            {
                if (ind[0] == 0)
                {

                    return distribution[0] * first_input_c1;
                }
                else if (ind[0] == 1)
                {
                    return distribution[0] * first_input_c2;
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
                    return distribution[0] * combined_output;
                }
                else
                {
                    std::cout << "Wrong Info: No such two check node, inds are: (" << ind[0] << ", " << ind[1] << "). plz check..." << std::endl;
                    return output_joint;
                }
            }
            else
            {
                std::cout << "Wrong Info: No such check node combinations, ind size is: " << ind.size() << ". plz check..." << std::endl;
                return output_joint;
            }

            //-----------------------------------------------------------------------------------------
        }
    }
    else
    {
        std::cout << "Neiter vari nor check, plz check again..." << std::endl;
        return output_joint;
    }
}

void ME2_PBRL_DE::type_distribution_update(std::vector<std::vector<double>> edge_distribution,
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
    if(strcmp(type,"check")==0)
    {
        llr_combined_output_joint = llr_combination(output_joint, 0.001);
        
    }
    else if(strcmp(type,"vari")==0)
    {
        llr_combined_output_joint = llr_combination(output_joint, 0.005);
    } 
    else if (strcmp(type,"check_min")==0)
    {
        llr_combined_output_joint = output_joint;
    }  
    clipped_cvd = clip_prob(llr_combined_output_joint, pow(10, -10.0));
    ave_joinprob_llr(clipped_cvd, pow(10.0, -80.0));
    llr = llr_cal(clipped_cvd);
    llr_file_name = std::string()+type+"_s" + std::to_string(socket) + "_llr_iteration_" + std::to_string(iter) + ".txt";
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
           // std::cout << "Info: Iteration " << iter << ", socket " << socket << ", variable mutual information update: " << IB_ins.mi << std::endl;
            std::cout<< iter<<"  " <<IB_ins.mi<<std::endl;
            break;
        default:
            std::cout << "Wrong Info: In function (type_distribution_update), type is vari, but no socket " << socket << ".. please check again." << std::endl;
            break;
        }
    }
    else if (strcmp(type, "check") == 0||strcmp(type, "check_min") == 0)
    {
        switch (socket)
        {
        case 0:
            check_representation_1.push_back(llr_cal(IB_ins.prob_join_xt));
            check_threshold_1.push_back(IB_ins.threshold);
            check_pmf_1 = IB_ins.prob_join_xt;
            //std::cout << "Info: Iteration " << iter << ", socket " << socket << ", check mutual information update: " << it_mi(IB_ins.prob_join_xt) << std::endl;            
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


bool ME2_PBRL_DE::read_decription()
{
    std::ifstream myfile(pbrl_met_description);
    if (myfile.is_open())
    {
        int cardi; //number of terms in each distribution

        //-------variable node -------------------------
        myfile >> cardi;
        variable_node_degree.resize(cardi, std::vector<double>(4, -1));
        for (unsigned ii = 0; ii < variable_node_degree.size(); ii++)
        {
            for (unsigned jj = 0; jj < 4; jj++)
            {
                myfile >> variable_node_degree[ii][jj];
            }
        }

        //------check node----------------------------
        myfile >> cardi;
        check_node_degree.resize(cardi, std::vector<double>(3, -1));
        for (unsigned ii = 0; ii < check_node_degree.size(); ii++)
        {
            for (unsigned jj = 0; jj < 3; jj++)
            {
                myfile >> check_node_degree[ii][jj];
            }
        }

        //------variable edge socket 1----------------------
        myfile >> cardi;
        vari_edge_deg_1.resize(cardi, std::vector<double>(4, -1));
        for (unsigned ii = 0; ii < vari_edge_deg_1.size(); ii++)
        {
            for (unsigned jj = 0; jj < 4; jj++)
            {
                myfile >> vari_edge_deg_1[ii][jj];
            }
        }

        //-------check edge socket 1---------------------------
        myfile >> cardi;
        check_edge_deg_1.resize(cardi, std::vector<double>(3, -1));
        for (unsigned ii = 0; ii < check_edge_deg_1.size(); ii++)
        {
            for (unsigned jj = 0; jj < 3; jj++)
            {
                myfile >> check_edge_deg_1[ii][jj];
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

int ME2_PBRL_DE::density_evolution()
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
    //-----check correctness of channel distirbution---------------------------
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
    //-------------------------------------------------------------------

    //-----push back channel threshold and reconstruction-------
    channel_threshold.push_back(channel_IB.threshold);  
    channel_representation.push_back(llr_cal(channel_IB.prob_join_xt));
    //-------------------------------------------------------

    vari_pmf_1 = np_channel_pmf;//channel_IB_1.prob_join_xt;
    vari_pmf_2 = np_channel_pmf;
    std::cout << "Info: Finished channel quantization ..." << std::endl;
    // std::vector<std::vector<double>> punc_check_edge_dist{{1.0/4,17,0},{1.0/4,6,1},{2.0/4,5,1}};
    // std::vector<std::vector<double>> punc_vari_edge_dist{{0.2667,-1,3,0},{0.7333,-1,4,0}};
    std::vector<std::vector<double>> punc_check_edge_dist=check_edge_deg_1;
    std::vector<std::vector<double>> punc_vari_edge_dist=vari_edge_deg_1;

    //---------------------density evolution------------------------------------
    for (unsigned ii = 0; ii < max_iter; ii++)
    {
        if (ii == 0)
        {
            //-----------------update check edge pmf---------------------------------
            type_distribution_update(punc_check_edge_dist, "check_min", ii, 0);
            //-----------------update variable node pmf------------------------------
            type_distribution_update(punc_vari_edge_dist, "vari", ii, 0);
        }
        else
        {
            //std::cout<<"---------------------------"<<std::endl;
            //-----------------update check edge pmf---------------------------------
            type_distribution_update(check_edge_deg_1, "check_min", ii, 0);
            //-----------------update variable node pmf------------------------------
            type_distribution_update(vari_edge_deg_1, "vari", ii, 0);
        }
    }
    RQF_output();
    return 1;
}

void ME2_PBRL_DE::RQF_output()
{
    //-------------------First Writeout channel_information---------------------
    std::string channelquan_file = "channel_quantizer_" + suffix + ".txt";
    std::string channelrec_file = "channel_reconstruction_" + suffix + ".txt";
    std::ofstream handle_channel_quantizer(channelquan_file);
    std::vector<std::vector<double>> cur_two_dim_vec;
    
    //------------------output channel information---------------------------------
    if (handle_channel_quantizer.is_open())
    {
        handle_channel_quantizer << quantization_size << "  " << std::endl;
        cur_two_dim_vec = channel_threshold;
        for (const auto &aa : cur_two_dim_vec)
        {
            for (const auto &bb : aa)
                handle_channel_quantizer << bb*sigma2/2.0 << "  ";
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
    threholdfile = "threshold_s1_" + suffix + ".txt";
    recons_file = "reconstruction_s1_" + suffix + ".txt";
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
        outpufile << std::endl;
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
        handle_reconstrction << std::endl;
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

}




 //-----punctured node probability distribution----------
    // temp_two_dim_vec.clear();
    // temp_two_dim_vec.push_back(channel_dist[0] * channel_observation[0]);
    // temp_two_dim_vec.push_back(channel_dist[0] * channel_observation[1]);
    // temp_two_dim_vec[0].push_back(channel_dist[1] * 0.5 * 1 / (1 + exp(small_offset)));
    // temp_two_dim_vec[0].push_back(channel_dist[1] * 0.5 * exp(small_offset) / (1 + exp(small_offset)));
    // temp_two_dim_vec[1].push_back(channel_dist[1] * 0.5 * exp(small_offset) / (1 + exp(small_offset)));
    // temp_two_dim_vec[1].push_back(channel_dist[1] * 0.5 * 1 / (1 + exp(small_offset)));
    // prob_sort(temp_two_dim_vec);
    // IB_kernel channel_IB_1(temp_two_dim_vec, quantization_size, ib_runtime);
    // channel_IB_1.Progressive_MMI();
    //-------------------------------------------------------

    //channel_threshold.push_back(channel_IB_1.threshold);
    // channel_representation.push_back(llr_cal(channel_IB_1.prob_join_xt));


        /*if (strcmp(type, "vari") == 0)
    {
        if (distribution.size() == 4)
        {
            if (distribution[1] > 0) //>0 not pucture
            {
                first_input = np_channel_pmf;
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
            }
            else //<0 puncture
            {
                if (distribution[2] > 0)
                {
                    first_input=llr_permutation(check_pmf_1,0.001);
                    second_input = check_pmf_1;
                    for (int ii = 0; ii < int(distribution[2])-1; ii++)
                    {
                        combined_output = prob_combination(first_input, second_input, type);
                        prob_sort(combined_output);
                        prob_llr_combined = llr_combination(combined_output, 0.001);
                        ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
                        first_input = prob_llr_combined;
                    }
                }
            }
            temp = distribution[0] * first_input[0];
            std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[0]));
            temp = distribution[0] * first_input[1];
            std::copy(temp.begin(), temp.end(), std::back_inserter(output_joint[1]));
            prob_sort(output_joint);
            return output_joint;
        }
        else
        {
            std::cout << "Wrong Info: It's not a variable node distribution vector (length: " << distribution.size() << "), plz check again" << std::endl;
            return output_joint;
        }
    }*/