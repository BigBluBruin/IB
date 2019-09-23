#pragma once
#include <vector>
#include <math.h>
#include <functional>
#include "Information_Bottleneck/overloadvec.h"
//This program gives a list of mutual inforamtion toolbox

double it_entropy(std::vector<double> dist);
double it_mi(std::vector<std::vector<double>> joindist);
double it_kl(std::vector<double> appdist, std::vector<double> trudist);
double it_js(double p1, std::vector<double> dist1, double p2, std::vector<double> dist2);
void ave_prob(std::vector <double> & prob);
void ave_joinprob(std::vector<std::vector<double>> & joinprob);
void flow_stop(double & in);