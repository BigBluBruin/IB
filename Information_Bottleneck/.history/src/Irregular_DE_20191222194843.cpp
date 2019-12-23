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
    std::vector<std::vector<double>> combined_check_dist(2), combined_vari_dist(2),llr_combined_check_dist, llr_combined_vari_dist, clipped_ccd,clipped_cvd;
    std::vector<double> temp;
    double most_left=-5.0;
    double most_right=5.0;
    int partition_number=3000;
    std::vector<std::vector<double>> channel_observation = gaussian_disretization(most_left, most_right, partition_number, sigma2);
    IB_kernel channel_IB(channel_observation, quantization_size, ib_runtime);
    channel_IB.smIB();
    std::vector<std::vector<double>> first_input = channel_IB.prob_join_xt;
    std::vector<std::vector<double>> second_input = first_input;
    std::vector<std::vector<double>> combined_output, prob_llr_combined;
    std::cout<<"finished channel quantization...."<<std::endl;

    //-------------------start iteration------------------------------------------
    for (unsigned iter = 0; iter < max_iter; iter++)
    {
        combined_check_dist[0].clear();
        combined_check_dist[1].clear();
        combined_vari_dist[0].clear();
        combined_vari_dist[1].clear();


        //---------------check node--------------------
        for (unsigned ii = 0; ii < check_edge_dist.size() - 2; ii++)
        {
            combined_output = prob_combination(first_input, second_input, "check");
            //combined_output = prob_combination_v2(first_input, second_input, "check",pow(10,-10));
            prob_sort(combined_output);
            prob_llr_combined = llr_combination(combined_output, 0.001);
            ave_joinprob_llr(prob_llr_combined, pow(10.0, -80.0));
            first_input = prob_llr_combined;
            
            if(check_edge_dist[ii+2]!=0)
            {
                temp=check_edge_dist[ii+2]*first_input[0];
                std::copy(temp.begin(),temp.end(),std::back_inserter(combined_check_dist[0]));
                std::cout<<"has been here"<<std::endl;
                temp=check_edge_dist[ii+2]*first_input[1];
                std::copy(temp.begin(),temp.end(),std::back_inserter(combined_check_dist[1]));
            }
        }
        prob_sort(combined_check_dist);
        llr_combined_check_dist = llr_combination(combined_check_dist, 0.001);
        clipped_ccd=clip_prob(llr_combined_check_dist,pow(10,-14.0));
        ave_joinprob_llr(clipped_ccd, pow(10.0, -80.0));
        std::vector<double> llr = llr_cal(clipped_ccd);
        std::string llr_file_name = "check_llr_iteration_" + std::to_string(iter) + ".txt";
        std::ofstream llr_file(llr_file_name);
        if (llr_file.is_open())
        {
            for (unsigned index = 0; index < llr.size(); index++)
            {
                llr_file << llr[index] << "  " << clipped_ccd[0][index] + clipped_ccd[1][index] << std::endl;
            }
        }
        llr_file.close();
        IB_kernel check_IB(clipped_ccd, quantization_size, ib_runtime);
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
                std::cout<<combined_vari_dist[1].size()<<"  "<<std::endl;
            }
           
        }
        prob_sort(combined_vari_dist);
        if(it_mi(combined_vari_dist)>0.99)
            llr_combination_interval=1;
        llr_combined_vari_dist = llr_combination(combined_vari_dist,llr_combination_interval);
        clipped_cvd=clip_prob(llr_combined_vari_dist,pow(10,-14.0));
        ave_joinprob_llr(clipped_cvd, pow(10.0, -80.0));
        llr = llr_cal(clipped_cvd);
        llr_file_name = "vari_llr_iteration_" + std::to_string(iter) + ".txt";
        llr_file.open(llr_file_name);
        if (llr_file.is_open())
        {
            for (unsigned index = 0; index < llr.size(); index++)
            {
                llr_file << llr[index] << "  " << clipped_cvd[0][index] + clipped_cvd[1][index] << std::endl;
            }
        }
        llr_file.close();
        IB_kernel vari_IB(clipped_cvd, quantization_size, ib_runtime);
        vari_IB.smIB();
        vari_representation.push_back(llr_cal(vari_IB.prob_join_xt));
        vari_threshold.push_back(vari_IB.threshold);
        //------------------------------------------------------------------------

        //---------------check node------------------------
        first_input = vari_IB.prob_join_xt;
        second_input = first_input;

       
        //--------------Output Info------------------------
        std::cout << "DE Info:  iter--" << iter + 1 << "--mi--" << vari_IB.mi <<"---true mi ---"<< it_mi(clipped_cvd)<<std::endl;
        std::cout << "Variable node threshold:" << std::endl;
        for (const auto &aa : vari_IB.threshold)
            std::cout << aa << "  ";
        std::cout << std::endl;
        std::cout << "cluster Info  :" << std::endl;
        for (const auto &aa : vari_IB.cluster)
            std::cout << aa << "  ";
        if (1 - vari_IB.mi < stop_threshold)
        {
            if (iter < max_iter - 3)
            {
                std::cout << "converge to fast..." << std::endl;
                //return 1;
            }
            else
            {
                if (iter == max_iter - 1)
                {
                    std::cout << "find threshold" << std::endl;
                    //return 3;
                }
            }
        }
        else
        {
            if (iter == max_iter - 1)
            {

                std::cout << "converge to slow..." << std::endl;
                //return 2;
            }
        }
    }


     //---------------file out info---------------------
        //Part 1: threholds 
        std::ofstream outpufile("threshold.txt");
        if (outpufile.is_open())
        {
            outpufile<<quantization_size<<" "<<std::endl;
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
            outpufile.close();
        }

        //Part 2: Reconstruction function
        std::ofstream handle_reconstrction("reconstruction.txt");
        if (handle_reconstrction.is_open())
        {
            handle_reconstrction << quantization_size << "  " << std::endl;
            for (const auto &iter : check_representation)
            {
                for (const auto &term : iter)
                {
                    handle_reconstrction << term << "  ";
                }
                handle_reconstrction << std::endl;
            }
            handle_reconstrction << std::endl
                                 << std::endl;
            for (const auto &iter : vari_representation)
            {
                for (const auto &term : iter)
                {
                    handle_reconstrction << term << "  ";
                }
                handle_reconstrction << std::endl;
            }
            handle_reconstrction.close();
        }
        else
        {
            std::cout << "Reconstruction failed to wirte ..." << std::endl;
        }

        //Part 3: channel quantizer
        std::ofstream handle_channel_quantizer("channel_quantizer.txt");
        if (handle_channel_quantizer.is_open())
        {
            int cur=0;
            handle_channel_quantizer << most_left << "  " << most_right << "  " << partition_number << std::endl;
            for (unsigned index = 0; index < channel_IB.cluster.size(); index++)
            {
                for (unsigned index2 = 0; index2 < channel_IB.cluster[index]; index2++)
                {
                    handle_channel_quantizer<<cur<<"  ";
                }
                cur=cur+1;
            }
            handle_channel_quantizer.close();
        }
        else
        {
            std::cout << "Fail to write channel quantizer ..." << std::endl;
        }

        //Part 4: channel Reconstruction
        std::ofstream handle_channel_recons("channel_reconstruction.txt");
        if (handle_channel_recons.is_open())
        {
            handle_channel_recons << quantization_size << "  ";
            for (unsigned index = 0; index < quantization_size; index++)
            {
                handle_channel_recons<<log(channel_IB.prob_join_xt[0][index]/channel_IB.prob_join_xt[1][index])<<"  "<<std::endl;
            }
            handle_channel_recons.close();
        }
        else
        {
            std::cout << "Fail to write channel reconstruction file ..." << std::endl;
        }

    return 1;

}