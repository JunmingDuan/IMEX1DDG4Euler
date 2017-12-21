#include "DGFEMSpace1D_GSL.h"

VEC<double> DGFEMSpace1D::LeftBoundary(const SOL& sol, const VEC<double>& U0,
    int* ghost_i, double* ghost_local_pnt, const double t) {
  VEC<double> U(DIM);
  if(BDL == 0) {
    U.resize(DIM, 0);
    *ghost_i = 0;
    *ghost_local_pnt = 1;
  }
  else if(BDL == 1) {
    U = U0;
    *ghost_i = 0;
    *ghost_local_pnt = -1;
  }
  else if(BDL == 2) {
    *ghost_i = Nx-1;
    *ghost_local_pnt = 1;
    U = Composition(sol,*ghost_i,mesh[Nx],t);
  }
  return U;
}

VEC<double> DGFEMSpace1D::RightBoundary(const SOL& sol, const VEC<double>& U0,
    int* ghost_i, double* ghost_local_pnt, const double t) {
  VEC<double> U(DIM);
  if(BDR == 0) {
    U.resize(DIM, 0);
    *ghost_i = Nx-1;
    *ghost_local_pnt = -1;
  }
  else if(BDR == 1) {
    U = U0;
    *ghost_i = Nx-1;
    *ghost_local_pnt = 1;
  }
  else if(BDR == 2) {
    *ghost_i = 0;
    *ghost_local_pnt = -1;
    U = Composition(sol,*ghost_i,mesh[0],t);
  }
  return U;
}

