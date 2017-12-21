#include "scaling_limiter.h"

void scaling_limiter::run(const VEC<VECD>& av, const VEC<VEC<VECD>>& val, SOL& sol)
{
  u_int n = av.size();
  double t, min_val;

  for(u_int i = 0; i < n; ++i) {
    u_int m = val[i].size();
    //limit rho
    min_val = EPS;
    for(u_int j = 0; j < m; ++j) {
      if(val[i][j][0] < min_val) min_val = val[i][j][0];
    }
    if(min_val < EPS) {
      std::cout << "Use PP limiter on rho!! i: " << i << std::endl;
      t = (av[i][0] - EPS) / (av[i][0] - min_val);
      //modify moments whos order > 1
      for(u_int k = 1; k < K; ++k) {
        sol[i][k][0] *= t;
      }
    }
  }
  if(DIM == 3) {
    //limit internal energy
    for(u_int i = 0; i < n; ++i) {
      u_int m = val[i].size();
      double ei = av[i][2] - 0.5*pow(av[i][1],2)/av[i][0];
      double val_e;
      min_val = EPS;
      for(u_int j = 0; j < m; ++j) {
        val_e = val[i][j][2] - 0.5*pow(val[i][j][1],2)/val[i][j][0];
        if(val_e < min_val) min_val = val_e;
      }
      if(min_val < EPS) {
        std::cout << "Use PP limiter on e!! i: " << i << std::endl;
        t = (ei - EPS) / (ei - min_val);
        //modify moments whos order > 1
        for(u_int k = 1; k < K; ++k) {
          sol[i][k] *= t;
        }
      }
    }
  }
}

