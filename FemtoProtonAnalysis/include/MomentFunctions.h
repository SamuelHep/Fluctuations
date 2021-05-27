#ifndef MOMENT_FUNCTIONS_H
#define MOMENT_FUNCTIONS_H

#include <vector>
#include <utility>

long double q_r_s(std::vector<std::pair<int,double>> &track_eff_pair, int r, int s);

std::vector<std::vector<long double>> make_all_q_s( std::vector<std::pair<int,double>> & track_eff_pair_vec, int rMax, int sMax );


#endif
