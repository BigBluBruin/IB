#include "Quantized_BP_PBRL.h"


Quantized_BP_PBRL::Quantized_BP_PBRL(std::string H_filename, std::string Pbrl_description_filename ,
                std::string Quantizer_name, std::vector<double> Parameters,
                int Target_error, std::string Channel_recon_file, std::string Recon_file,
                std::string Threshold_file)
{
    h_filename=H_filename;
    pbrl_description_filename = Pbrl_description_filename;
    cq_ins.quantizer_file = Quantizer_name;
    parameters = Parameters;
    target_error = Target_error;
    total_frames.assign(parameters.size(), -1);
    total_iterations.assign(parameters.size(), -1);
    channel_recon_file = Channel_recon_file;
    recon_file = Recon_file;
    threshold_file = Threshold_file;
}

void Quantized_BP_PBRL::fill_in()
{
    //---------------parity check matrix---------------------------------
    h_ins.filename = h_filename;
    h_ins.pbrl_description_file = pbrl_description_filename;
    h_ins.Read_Parity_Check_Matrix();
    h_ins.read_pbrl_info();
    h_ins.assign_socket();

    //---------quantization info----------------------------
    cq_ins.read_channel_quantizer();
    //-------------------------------------------------------

    //---------channel reconstruction------------------------
    std::ifstream handle_channel_rec(channel_recon_file);
    unsigned int quan_size;
    if (handle_channel_rec.is_open())
    {
        handle_channel_rec >> quan_size;
        channel_recons.assign(quan_size, -1.0);
        for (unsigned index = 0; index < quan_size; index++)
            handle_channel_rec >> channel_recons[index];
        std::cout << "Channel reconstruction file Successfull..."
                  << std::endl;
    }
    else
    {
        std::cout << "Channel reconstruction file NOT found"
                  << std::endl;
    }    
    //-------------------------------------------------------

    //--------reconstruction file-----------------------------
    std::ifstream handel_reconstruction(recon_file);
    if (handel_reconstruction.is_open())
    {
        handel_reconstruction >> quan_size >> max_iter;
        check_recons.resize(max_iter);
        vari_recons.resize(max_iter);
        for (int index1 = 0; index1 < max_iter; index1++)
        {
            check_recons[index1].assign(quan_size, -1);
            for (unsigned index2 = 0; index2 < quan_size; index2++)
            {
                handel_reconstruction >> check_recons[index1][index2];
            }
        }
        for (int index1 = 0; index1 < max_iter; index1++)
        {
            vari_recons[index1].assign(quan_size, -1);
            for (unsigned index2 = 0; index2 < quan_size; index2++)
            {
                handel_reconstruction >> vari_recons[index1][index2];
            }
        }
        std::cout << "Reconstruction file Successfull..."
                  << std::endl;

    }
    else
    {
        std::cout << "Reconstruction file NOT found"
                  << std::endl;
    }

    //--------------------------------------------------------

    //--------threshold files---------------------------------
    std::ifstream handle_quan_threshold(threshold_file);
    if (handle_quan_threshold.is_open())
    {
        handle_quan_threshold >> quan_size >> max_iter;
        check_threshold.assign(max_iter, std::vector<double>(1, -1));
        vari_threshold.assign(max_iter, std::vector<double>(1, -1));
        for (int index1 = 0; index1 < max_iter; index1++)
        {
            check_threshold[index1].assign(quan_size / 2 - 1, -1);
            for (unsigned index2 = 0; index2 < quan_size / 2 - 1; index2++)
            {
                handle_quan_threshold >> check_threshold[index1][index2];
            }
        }
        for (int index1 = 0; index1 < max_iter; index1++)
        {
            vari_threshold[index1].assign(quan_size / 2 - 1, -1);
            for (unsigned index2 = 0; index2 < quan_size / 2 - 1; index2++)
            {
                handle_quan_threshold >> vari_threshold[index1][index2];
            }
        }
        std::cout << "Threshold file Successfull..."
                  << std::endl;
    }
    else
    {
        std::cout << "Threshold file NOT found"
                  << std::endl;
    }
    //--------------------------------------------------------

}

void Quantized_BP_PBRL::noise_generator(std::vector<double> &cwds, double parameter)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> distribution(0, sqrt(parameter));
    for (unsigned inc = 0; inc < cwds.size(); inc++)
    {
        cwds[inc] += distribution(generator);
    }
}

unsigned Quantized_BP_PBRL::threshold_quantization(std::vector<double> &threshold, double value)
{
    unsigned ind = 10000;
    unsigned quan_size = threshold.size() + 1;
    bool negtive = true;
    if (value > 0)
    {
        negtive = false;
        value = -1 * value;
    }
    if (value <= threshold[0])
    {
        ind = 0;
    }
    else if (value > threshold[quan_size - 2])
    {
        ind = quan_size - 1;
    }
    else
    {
        for (unsigned index = 0; index < quan_size - 2; index++)
        {
            if (value > threshold[index] && value <= threshold[index + 1])
            {
                ind = index + 1;
            }
        }
    }
    if (!negtive)
    {
        ind = 2 * quan_size - 1 - ind;
    }
    return ind;
}

bool Quantized_BP_PBRL::iscwds(std::vector<int> final_bits)
{
    std::vector<int> edge_bits(h_ins.edge_num, -1);
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            edge_bits[h_ins.edge_v[ii][jj]] = final_bits[ii];
        }
    }

    for (int ii = 0; ii < h_ins.check_thre; ii++)
    {
        int checks = 0;
        for (int jj = 0; jj < h_ins.check_degreetable[ii]; jj++)
        {
            checks = (checks + edge_bits[h_ins.edge_c[ii][jj]]) % 2;
        }
        if (checks > 0)
        {
            return false;
        }
    }

    return true;
}

bool Quantized_BP_PBRL::decoder(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits)
{
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<double> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    unsigned quantized_val;
    double full_return;
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    double final_codewords;
    std::vector<double> message;
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        if(ii<h_ins.punc_num)
        {
            rx[ii] = -1; 
        }
        else
        {
            rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
        }
        
        
    }
    //Ini v2c messages
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            msg_v2c[h_ins.edge_v[ii][jj]] = rx[ii];
        }
    }

    //Start Iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        //c2v update
        for (int ii = 0; ii < h_ins.check_num; ii++)
        {
            int cur_dc = h_ins.check_degreetable[ii];
            for (int jj = 0; jj < cur_dc; jj++)
            {
                message.clear();
                for (int kk = 0; kk < cur_dc; kk++)
                {
                    if (kk != jj && msg_v2c[h_ins.edge_c[ii][kk]]!= -1)
                    {
                        message.push_back(msg_v2c[h_ins.edge_c[ii][kk]]);
                    }
                }
                if (message.size() != cur_dc - 1)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]] = -1;
                }
                else
                {
                    full_return = check_node_operation(message);
                    quantized_val = threshold_quantization(check_threshold[cur_iter], full_return);
                    if (quantized_val != 10000)
                    {
                        msg_c2v[h_ins.edge_c[ii][jj]] = check_recons[cur_iter][quantized_val];
                        //std::cout<<msg_c2v[h_ins.edge_c[ii][jj]]<<std::endl;
                    }
                    else
                    {
                        std::cout << "variable node quantization failed..." << std::endl;
                    }
                }
            }
        }

        //v2c update
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii]; //find current dv
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message.clear();
                if(ii>=h_ins.punc_num)
                {
                    message.push_back(rx[ii]);
                }
                //std::cout<<"have been 2.5"<<std::endl;
                for (int kk = 0; kk < cur_dv; kk++)
                {
                    //collect data
                    if (kk != jj && msg_c2v[h_ins.edge_v[ii][kk]] != -1)
                    {
                        message.push_back(msg_c2v[h_ins.edge_v[ii][kk]]);
                    }
                }
                full_return = vari_node_operation(message);
                quantized_val = threshold_quantization(vari_threshold[cur_iter], full_return);
                if (quantized_val != 10000)
                {
                    msg_v2c[h_ins.edge_v[ii][jj]] = vari_recons[cur_iter][quantized_val];
                }
                else
                {
                    std::cout << "variable node quantization failed..." << std::endl;
                }
            }
        }

        //final decision
        for (int ii = 0; ii < h_ins.vari_thre; ii++)
        {

            int cur_dv = h_ins.vari_degreetable[ii];
            message.clear();
            if (ii >= h_ins.punc_num)
            {
                message.push_back(rx[ii]);
            }
            for (int jj = 0; jj < cur_dv; jj++)
            {
                if(msg_c2v[h_ins.edge_v[ii][jj]] != -1)
                {
                     message.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
                }
            }
            final_codewords = vari_node_operation(message);
            if (final_codewords > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }

        for (int ii = h_ins.vari_thre; ii < h_ins.vari_num; ii++)
        {
            if (rx[ii] > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }

        //check sum
        if (iscwds(final_bits))
        {
            iteration = iteration + cur_iter + 1;
            int sum = 0;
            for (int ii = 0; ii < h_ins.vari_thre; ii++)
            {
                sum = sum + final_bits[ii];
            }
            if (sum == 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    iteration = iteration + max_iter;
    return false;
}

void Quantized_BP_PBRL::main_simulation(const char ind[],const char suffix[])
{
    int total_num = 0;
    int error_num = 0;
    int iteration;
    double cur_para;
    std::vector<double> codewords;
    std::vector<int> final_bits;
    for (unsigned ii = 0; ii < parameters.size(); ii++)
    {
        iteration = 0;
        cur_para = pow(10, (-0.1 * parameters[ii]) / (2.0 * h_ins.rate));
        std::string result_file = "Result_PBRL_n_" + std::to_string(h_ins.vari_num) + "_k_" + std::to_string(h_ins.check_num) + "_den_" + suffix + "_Para_" + std::to_string(parameters[ii]) + "_ind_" + ind + ".txt";
        std::ofstream result_bar;
        do
        {
            total_num = total_num + 1;
            codewords.assign(h_ins.vari_num, 1);
            final_bits.assign(h_ins.vari_num, -1);
            //add noise
            noise_generator(codewords, cur_para);
            bool result = decoder(codewords, iteration, final_bits);
            if (!result)
            {
                /* fail */
                error_num = error_num + 1;
                std::cout << "Para: " << parameters[ii] << " error: " << error_num << "fer: " << (double)error_num / (double)total_num << std::endl;
                //you can also output failed codewords here....
            }
            if (total_num % 100 == 1)
            {
                result_bar.open(result_file);
                result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
                result_bar.close();
            }

        } while (error_num < target_error);
        result_bar.open(result_file);
        result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
        result_bar.close();
    }
}