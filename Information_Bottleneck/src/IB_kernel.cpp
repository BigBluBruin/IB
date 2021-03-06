#include "Information_Bottleneck/IB_kernel.h"

IB_kernel::IB_kernel(std::vector<std::vector<double>> input, unsigned quan, int Max_run)
{
    prob_join_xy = input;
    quan_size = quan;
    max_run = Max_run;
    prob_join_xt.assign(2,std::vector<double>(quan_size,0));
    prob_t.assign(quan_size,0);
    threshold.assign(quan_size/2-1,-1);
}

std::vector<std::vector<double>> IB_kernel::quantize_to_xt(std::vector<std::vector<double>> &input, std::vector<unsigned> &cluster)
{
    std::vector<std::vector<double>> join_prob_xt(2, std::vector<double>(cluster.size(), 0.0));
    std::vector<double>::iterator iter1 = input[0].begin();
    std::vector<double>::iterator iter2 = input[1].begin();
    int i = 0;
    for (const unsigned &cluval : cluster)
    {
        join_prob_xt[0][i] = std::accumulate(iter1, iter1 + cluval, 0.0);
        join_prob_xt[1][i] = std::accumulate(iter2, iter2 + cluval, 0.0);
        iter1 = iter1 + cluval;
        iter2 = iter2 + cluval;
        i=i+1;
    }
    return join_prob_xt;
}

std::vector<unsigned> IB_kernel::random_cluster(const unsigned total_num, const unsigned quan_size)
{

    //initialize random seed
    std::random_device rd;
    std::mt19937 generator(rd());
    //srand(time(0));

    //Para Initial
    std::vector<unsigned> remaining(total_num - 1);
    for (unsigned ii = 0; ii < total_num - 1; ii++)
    {
        remaining[ii] = ii + 1;
    }


    std::vector<unsigned> chosen(quan_size - 1);
    int selected;
    //start
    for (unsigned ii = 0; ii < quan_size - 1; ii++)
    {
        std::uniform_int_distribution<> distribution(0, total_num - 2 - ii);
        selected = distribution(generator);
        chosen[ii] = remaining[selected];
        remaining.erase(remaining.begin() + selected);
    }
    chosen.push_back(0);
    chosen.push_back(total_num);
    std::sort(chosen.begin(), chosen.end());

    std::vector<unsigned> cluster(quan_size);
    for (unsigned ii = 0; ii < quan_size; ii++)
    {
        cluster[ii] = chosen[ii + 1] - chosen[ii];
    }


    return cluster;
}

void IB_kernel::smIB()
{
    //initialization...
    std::vector<double> left_join, right_join, sin_join;
    double left_cost, right_cost;
    double py, pt;
    bool finish = false;
    double best_mi=0;
    std::vector<unsigned> best_partition;
    std::vector<std::vector<double>> temp_join_xt;
    for (int run_ind = 0; run_ind < max_run; run_ind++)
    {
        //step 1: generate a random cluster
        std::vector<unsigned> left_par = random_cluster(prob_join_xy[1].size()/2, quan_size / 2);
        std::vector<double>::iterator iter1, iter2;
        unsigned counter;
        std::vector<unsigned> right_par = left_par;
        std::reverse(right_par.begin(), right_par.end());
        std::vector<unsigned> partition;
        std::copy(left_par.begin(), left_par.end(), std::back_inserter(partition));
        std::copy(right_par.begin(), right_par.end(), std::back_inserter(partition));
        bool left_done=false;
        bool right_done=false;
        /*back_inserter is a convenience function template that constructs a 
        std::back_insert_iterator for the container c with the type deduced from 
        the type of the argument.*/
        //------debugging-------------------
        /*for (const auto& term: partition)
            std::cout<<term<<std::endl;
        std::cout<<std::endl;*/
        //-----------------------------------
        do
        {
            left_done = true;
            right_done = true;
            //step 2: for all left part
            for (unsigned ii = 0; ii < quan_size / 2 - 1; ii++)
            {

                do
                { //each cluater has at least one element !
                    if (partition[ii] > 1)
                    {

                        //get left, right and single probabilities
                        /*x += y is equivalent to x = x + y*/
                        counter = 0;
                        for (unsigned jj = 0; jj < ii; jj++)
                            counter += partition[jj];
                        //left
                        iter1 = prob_join_xy[0].begin() + counter;
                        iter2 = iter1 + partition[ii] - 1;
                        left_join.resize(2);
                        left_join[0] = std::accumulate(iter1, iter2, 0.0);
                        iter1 = prob_join_xy[1].begin() + counter;
                        iter2 = iter1 + partition[ii] - 1;
                        left_join[1] = std::accumulate(iter1, iter2, 0.0);
                        //single
                        sin_join = {prob_join_xy[0][counter + partition[ii] - 1], prob_join_xy[1][counter + partition[ii] - 1]};
                        //right
                        iter1 = prob_join_xy[0].begin() + counter + partition[ii];
                        iter2 = iter1 + partition[ii + 1];
                        right_join.resize(2);
                        right_join[0] = std::accumulate(iter1, iter2, 0.0);
                        iter1 = prob_join_xy[1].begin() + counter + partition[ii];
                        iter2 = iter1 + partition[ii + 1];
                        right_join[1] = std::accumulate(iter1, iter2, 0.0);
                        py = left_join[0] + left_join[1];
                        pt = sin_join[0] + sin_join[1];
                        ave_prob(left_join);
                        ave_prob(sin_join);
                        //std::cout<<"left: "<<py<<"single: "<<pt<<" ";
                        left_cost = (py + pt) * it_js(py / (py + pt), left_join, pt / (py + pt), sin_join);
                        py = right_join[0] + right_join[1];
                        ave_prob(right_join);
                        right_cost = (py + pt) * it_js(py / (py + pt), right_join, pt / (py + pt), sin_join);
                        //std::cout<<"right: "<<py<<"single: "<<pt<<" ";
                        if (left_cost < right_cost)
                            finish = true;
                        else
                        {
                            left_done = false;
                            partition[ii]--;
                            partition[ii + 1]++;
                            partition[quan_size - 1 - ii]--;
                            partition[quan_size - 2 - ii]++;
                            //std::cout<<ii<<"th size is: "<<partition[ii]<<std::endl;
                        }
                        //std::cout<<left_cost<<"  "<<right_cost<<std::endl;
                    }
                    else
                    {
                        finish = true;
                    }

                } while (!finish);
            }
            //std::cout<<"left finished"<<std::endl;

            finish = false;
            //step 3: for all right part
            for (unsigned ii = 0; ii < quan_size / 2 - 1; ii++)
            {
                do
                { //each cluater has at least one element !
                    if (partition[ii + 1] > 1)
                    {

                        //get left, right and single probabilities
                        counter = 0;
                        for (unsigned jj = 0; jj < ii; jj++)
                            counter += partition[jj];
                        //left
                        iter1 = prob_join_xy[0].begin() + counter;
                        iter2 = iter1 + partition[ii];
                        left_join.resize(2);
                        left_join[0] = std::accumulate(iter1, iter2, 0.0);
                        iter1 = prob_join_xy[1].begin() + counter;
                        iter2 = iter1 + partition[ii];
                        left_join[1] = std::accumulate(iter1, iter2, 0.0);
                        //std::cout<<"finished left"<<std::endl;
                        //single
                        sin_join = {prob_join_xy[0][counter + partition[ii]], prob_join_xy[1][counter + partition[ii]]};
                        //right
                        iter1 = prob_join_xy[0].begin() + counter + partition[ii] + 1;
                        iter2 = prob_join_xy[0].begin() + counter + partition[ii] + partition[ii + 1];
                        right_join.resize(2);
                        right_join[0] = std::accumulate(iter1, iter2, 0.0);
                        iter1 = prob_join_xy[1].begin() + counter + partition[ii] + 1;
                        iter2 = prob_join_xy[1].begin() + counter + partition[ii] + partition[ii + 1];
                        right_join[1] = std::accumulate(iter1, iter2, 0.0);
                        //std::cout<<"finished right"<<std::endl;
                        py = left_join[0] + left_join[1];
                        pt = sin_join[0] + sin_join[1];
                        //calculate
                        ave_prob(left_join);
                        ave_prob(sin_join);
                        left_cost = (py + pt) * it_js(py / (py + pt), left_join, pt / (py + pt), sin_join);
                        py = right_join[0] + right_join[1];
                        ave_prob(right_join);
                        right_cost = (py + pt) * it_js(py / (py + pt), right_join, pt / (py + pt), sin_join);
                        //std::cout<<left_cost<<"  "<<right_cost<<"size: "<<std::endl;
                        if (left_cost > right_cost)
                            finish = true;
                        else
                        {
                            right_done = false;
                            partition[ii]++;
                            partition[ii + 1]--;
                            partition[quan_size - 1 - ii]++;
                            partition[quan_size - 2 - ii]--;
                        }
                    }
                    else
                    {
                        finish = true;
                    }

                } while (!finish);
            }
        } while (!(left_done && right_done));

        //std::cout<<"finish left and right in "<<run_ind<<std::endl;
        //step 4: store and replace
        if(run_ind==0)
        {
            temp_join_xt=quantize_to_xt(prob_join_xy,partition);
            best_mi=it_mi(temp_join_xt);
            best_partition=partition;
        }
        else
        {

            temp_join_xt=quantize_to_xt(prob_join_xy,partition);
            if(it_mi(temp_join_xt)>best_mi)
            {
                best_mi=it_mi(temp_join_xt);
                best_partition=partition;
            }
            //std::cout<<it_mi(temp_join_xt)<<std::endl;
            if(isnan(it_mi(temp_join_xt)))
            {
                std::cout<<"------------------"<<std::endl;
                for(const auto & layer1: temp_join_xt)
                {
                    for(const auto& term: layer1)
                    {
                        std::cout<<term<<"  ";
                    }
                    std::cout<<std::endl;
                    
                }
                std::cout<<"******************"<<std::endl;
                for(const auto & term: partition)
                    std::cout<<term<<",";
                std::cout<<std::endl;
                return;
            }

        }
    }

    cluster=best_partition;
    prob_join_xt=quantize_to_xt(prob_join_xy,cluster);
    mi=it_mi(prob_join_xt);
    unsigned counter;
    for(unsigned index=1;index<=quan_size/2-1;index++)
    {
        counter=std::accumulate(cluster.begin(),cluster.begin()+index,0.0);
        threshold[index-1]=log(prob_join_xy[0][counter]/prob_join_xy[1][counter]);
    }
    for(unsigned ind=0;ind<quan_size;ind++)
    {
        prob_t[ind]=prob_join_xt[0][ind]+prob_join_xt[1][ind];
    }
}





unsigned IB_kernel::find_threshold(unsigned left_most, unsigned right_most)
{
    // In this instance, we shall obey C++ rule such that left point is inclued and right
    //point is not included. This is beneficial for writting program

    //Initialization
    //std::cout<<"left most point--"<<left_most<<"--right most--"<<right_most<<std::endl;
    std::vector<double> right_join_prob(2), left_join_prob(2), single_join_prob(2);
    unsigned total_time = right_most - left_most - 2;
    unsigned single_point_pos;
    single_point_pos=left_most + 1;
    double py, pt;
    double left_cost,right_cost;
    
    for (unsigned index = 0; index < total_time; index++)
    {
        single_point_pos = left_most + 1 + index;
        left_join_prob[0] = std::accumulate(prob_join_xy[0].begin() + left_most, prob_join_xy[0].begin() + single_point_pos, 0.0);
        left_join_prob[1] = std::accumulate(prob_join_xy[1].begin() + left_most, prob_join_xy[1].begin() + single_point_pos, 0.0);
        right_join_prob[0] = std::accumulate(prob_join_xy[0].begin() + single_point_pos + 1, prob_join_xy[0].begin() + right_most, 0.0);
        right_join_prob[1] = std::accumulate(prob_join_xy[1].begin() + single_point_pos + 1, prob_join_xy[1].begin() + right_most, 0.0);
        single_join_prob[0] = prob_join_xy[0][single_point_pos];
        single_join_prob[1] = prob_join_xy[1][single_point_pos];
        py = left_join_prob[0] + left_join_prob[1];
        pt = single_join_prob[0] + single_join_prob[1];
        //std::cout<<py<<"  "<<pt<<std::endl;
        ave_prob(left_join_prob);
        ave_prob(single_join_prob);
        left_cost = (py + pt) * it_js(py / (py + pt), left_join_prob, pt / (py + pt), single_join_prob);
        py = right_join_prob[0] + right_join_prob[1];
        ave_prob(right_join_prob);
        right_cost = (py + pt) * it_js(py / (py + pt), right_join_prob, pt / (py + pt), single_join_prob);
        if (right_cost<left_cost)
        {
            return single_point_pos;
        }
    }
    return single_point_pos;
}

void IB_kernel::Progressive_MMI()
{
    unsigned total_size = prob_join_xy[1].size();
    unsigned sym_size = total_size / 2;
    unsigned outer_time = log2(quan_size) - 1;
    unsigned itercount = outer_time-1;
    int cur_left_most;
    int cur_right_most;
    unsigned inner_length;
    int m;
    std::vector<unsigned> partition,partition_reverse;
    std::vector<unsigned> boundary;
    boundary.push_back(0);
    boundary.push_back(sym_size);
    for (unsigned index = 0; index < outer_time; index++)
    {
        inner_length = boundary.size();
        for (unsigned index2 = 0; index2 < inner_length - 1; index2++)
        {
            if(itercount >=0)
            {
                m = pow(2,itercount);
            }
            else
            {
                std::cout<<"Wrong info: variable itercount should be non-negtive"<<std::endl;
            }
            
            cur_left_most = boundary[index2];
            cur_right_most = boundary[index2 + 1];
            boundary.push_back(find_threshold(cur_left_most+(m-1), cur_right_most-(m-1)));
        }
        std::sort(boundary.begin(), boundary.end());
        // std::cout<<m<<"|| ";
        // display(boundary);
        itercount--;     
    }
    
    for (unsigned index = 0; index < quan_size / 2; index++)
    {
        partition.push_back(boundary[index + 1] - boundary[index]);
    }
    partition_reverse = partition;
    std::reverse(partition_reverse.begin(), partition_reverse.end());
    std::copy(partition_reverse.begin(), partition_reverse.end(), std::back_inserter(partition));
    cluster=partition;
    //display(cluster);
    //-----------------------------------------------
    prob_join_xt = quantize_to_xt(prob_join_xy, cluster);
    mi = it_mi(prob_join_xt);
    unsigned counter;
    for (unsigned index = 1; index <= quan_size / 2 - 1; index++)
    {
        counter = std::accumulate(cluster.begin(), cluster.begin() + index, 0.0);
        threshold[index - 1] = log(prob_join_xy[0][counter] / prob_join_xy[1][counter]);
    }
    for (unsigned ind = 0; ind < quan_size; ind++)
    {
        prob_t[ind] = prob_join_xt[0][ind] + prob_join_xt[1][ind];
    }
}


void IB_kernel::external_force(std::vector<unsigned> Cluster)
{
    cluster = Cluster;
    prob_join_xt = quantize_to_xt(prob_join_xy, cluster);
    mi = it_mi(prob_join_xt);
    unsigned counter;
    for (unsigned index = 1; index <= quan_size / 2 - 1; index++)
    {
        counter = std::accumulate(cluster.begin(), cluster.begin() + index, 0.0);
        threshold[index - 1] = log(prob_join_xy[0][counter] / prob_join_xy[1][counter]);
    }
    for (unsigned ind = 0; ind < quan_size; ind++)
    {
        prob_t[ind] = prob_join_xt[0][ind] + prob_join_xt[1][ind];
    }
}