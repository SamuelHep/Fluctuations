#ifndef MOMENT_FUNCTIONS_H
#define MOMENT_FUNCTIONS_H

#include <vector>
#include <utility>

double q_r_s(std::vector<std::pair<int,double>> &track_eff_pair, int r, int s);

std::vector<std::vector<double>> make_all_q_s( std::vector<std::pair<int,double>> & track_eff_pair_vec, int rMax, int sMax );


#endif
