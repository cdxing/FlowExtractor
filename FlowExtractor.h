int Bin_Centrality_01 = 4;
int Bin_rap = 8;
const Double_t _sigmaRange = 5.; // Sigma of the Fitting range
const Double_t _y_CM = -2.02;
const Double_t _n_jkk = 1;

const int cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
const int cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
const double rap_low_phi[8] = {-1.0, -0.6, -0.3, -0.1, 0., 0.1, 0.3, 0.6};
const double rap_up_phi[8]  = {-0.6, -0.3, -0.1, 0.,  0.1, 0.3, 0.6, 1.0};
const Double_t a_d_int_range[4] ={
  0.99,
  1.014,
  1.026,
  1.09
};
