#include <cstring>
#include <iostream>
#include <utility>
#include <vector>
#include <math.h>

using namespace std;

long double q_r_s(vector<pair<int,double>> &track_eff_pair, int r, int s)
{
  long double q = 0;
  
  for ( auto &p : track_eff_pair )
    {
      if( p.second <= 0 ) continue;
      q += pow(p.first,r)/pow(p.second,s);
    }

  return q;
}

vector<vector<long double>> make_all_q_s( vector<pair<int,double>> & track_eff_pair_vec, int rMax, int sMax )
{
  //Init an empty vector, the first entry isnt used to make it the same as the factorial moment indicies
  vector<vector<long double>> q_r_s_vec(rMax+1, vector<long double>(sMax+1, 0.0));

  for (int r=1;r<=rMax;r++)
    {
      for (int s=1;s<=sMax;s++)
	{
	  q_r_s_vec[r][s] = q_r_s(track_eff_pair_vec,r,s);
	}
    }

  return q_r_s_vec;
  
}

/*
charge efficiency 


q[1][1] ... q[6][6]

*/
