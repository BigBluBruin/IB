#include "Information_Bottleneck/IB_kernel.h"

IB_kernel::IB_kernel(std::vector<std::vector<double>> input, unsigned quan, int Max_run)
{
    prob_join_xy = input;
    quan_size = quan;
    max_run = Max_run;
    prob_join_xt.assign(2,std::vector<double>(quan_size,0));
}
std::vector<std::vector<double>> IB_kernel::quantize_to_xt(std::vector<std::vector<double>> &input, std::vector<unsigned> &cluster)
{
    std::vector<std::vector<double>> join_prob_xt(2, std::vector<double>(cluster.size(), 0));
    std::vector<double>::iterator iter1 = input[0].begin();
    std::vector<double>::iterator iter2 = input[1].begin();
    int i = 0;
    for (const unsigned &cluval : cluster)
    {
        join_prob_xt[0][i] = std::accumulate(iter1, iter1 + cluval, 0);
        join_prob_xt[1][i] = std::accumulate(iter2, iter2 + cluval, 0);
        iter1 = iter1 + cluval;
        iter2 = iter2 + cluval;
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

    for (const auto &term : remaining)
        std::cout << term;
    std::cout << std::endl;
    std::vector<unsigned> chosen(quan_size - 1);
    int selected;

    //start
    for (unsigned ii = 0; ii < quan_size - 1; ii++)
    {
        std::uniform_int_distribution<> distribution(0, total_num - 2 - ii);
        selected = distribution(generator);
        std::cout << "round: " << ii << "select:" << selected;
        chosen[ii] = remaining[selected];
        std::cout << "number" << chosen[ii] << std::endl;
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
    double best_mi;
    std::vector<unsigned> best_partition;
    std::vector<std::vector<double>> temp_join_xt;
    for (int run_ind = 0; run_ind < max_run; run_ind++)
    {
        //step 1: generate a random cluster
        std::vector<unsigned> left_par = random_cluster(prob_join_xy[1].size(), quan_size / 2);
        std::vector<double>::iterator iter1, iter2;
        unsigned counter;
        std::vector<unsigned> right_par = left_par;
        std::reverse(right_par.begin(), right_par.end());
        std::vector<unsigned> partition;
        std::copy(left_par.begin(), left_par.end(), std::back_inserter(partition));
        std::copy(right_par.begin(), right_par.end(), std::back_inserter(partition));

        //step 2: for all left part
        for (unsigned ii = 0; ii < quan_size / 2 - 1; ii++)
        {

            do
            { //each cluater has at least one element !
                if (partition[ii] > 1)
                {

                    //get left, right and single probabilities
                    counter = 0;
                    for (unsigned jj = 0; jj < counter; jj++)
                        counter += partition[jj];
                    //left
                    iter1 = prob_join_xy[0].begin() + counter;
                    iter2 = iter1 + partition[ii] - 1;
                    left_join.resize(2);
                    left_join[0] = std::accumulate(iter1, iter2, 0);
                    iter1 = prob_join_xy[1].begin() + counter;
                    iter2 = iter1 + partition[ii] - 1;
                    left_join[1] = std::accumulate(iter1, iter2, 0);
                    //single
                    sin_join = {prob_join_xy[0][counter + partition[ii] - 1], prob_join_xy[1][counter + partition[ii] - 1]};
                    //right
                    iter1 = prob_join_xy[0].begin() + counter + partition[ii];
                    iter2 = iter1 + partition[ii + 1];
                    right_join.resize(2);
                    right_join[0] = std::accumulate(iter1, iter2, 0);
                    iter1 = prob_join_xy[1].begin() + counter + partition[ii];
                    iter2 = iter1 + partition[ii + 1];
                    right_join[1] = std::accumulate(iter1, iter2, 0);
                    py = left_join[0] + left_join[1];
                    pt = sin_join[0] + sin_join[1];
                    left_cost = (py + pt) * it_js(py/(py+pt), left_join, pt/(py+pt), sin_join);
                    py = right_join[0] + right_join[1];
                    right_cost = (py + pt) * it_js(py/(py+pt), right_join, pt/(py+pt), sin_join);
                    if (left_cost < right_cost)
                        finish = true;
                    else
                    {
                        partition[ii]--;
                        partition[ii + 1]++;
                        partition[quan_size-1-ii]--;
                        partition[quan_size-2-ii]++;
                    }
                }
                else
                {
                    finish = true;
                }

            } while (finish);
        }

        //step 3: for all right part
        for (unsigned ii = 0; ii < quan_size / 2 - 1; ii++)
        {

            do
            { //each cluater has at least one element !
                if (partition[ii+1] > 1)
                {

                    //get left, right and single probabilities
                    counter = 0;
                    for (unsigned jj = 0; jj < counter; jj++)
                        counter += partition[jj];
                    //left
                    iter1 = prob_join_xy[0].begin() + counter;
                    iter2 = iter1 + partition[ii];
                    left_join.resize(2);
                    left_join[0] = std::accumulate(iter1, iter2, 0);
                    iter1 = prob_join_xy[1].begin() + counter;
                    iter2 = iter1 + partition[ii];
                    left_join[1] = std::accumulate(iter1, iter2, 0);
                    //single
                    sin_join = {prob_join_xy[0][counter + partition[ii]], prob_join_xy[1][counter + partition[ii]]};
                    //right
                    iter1 = prob_join_xy[0].begin() + counter + partition[ii]+1;
                    iter2 = prob_join_xy[0].begin() + counter + partition[ii + 1];
                    right_join.resize(2);
                    right_join[0] = std::accumulate(iter1, iter2, 0);
                    iter1 = prob_join_xy[1].begin() + counter + partition[ii]+1;
                    iter2 = prob_join_xy[1].begin() + counter + partition[ii+1];
                    right_join[1] = std::accumulate(iter1, iter2, 0);
                    py = left_join[0] + left_join[1];
                    pt = sin_join[0] + sin_join[1];
                    left_cost = (py + pt) * it_js(py, left_join, pt, sin_join);
                    py = right_join[0] + right_join[1];
                    right_cost = (py + pt) * it_js(py, right_join, pt, sin_join);
                    if (left_cost > right_cost)
                        finish = true;
                    else
                    {
                        partition[ii]++;
                        partition[ii + 1]--;
                        partition[quan_size-1-ii]++;
                        partition[quan_size-2-ii]--;
                    }
                }
                else
                {
                    finish = true;
                }

            } while (finish);
        }
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
        }
    }
    cluster=best_partition;
    prob_join_xt=quantize_to_xt(prob_join_xy,cluster);
    mi=it_mi(prob_join_xt);
    for(unsigned ind=0;ind<quan_size;ind++)
    {
        prob_t[ind]=prob_join_xt[0][ind]+prob_join_xt[1][ind];
    }
}

