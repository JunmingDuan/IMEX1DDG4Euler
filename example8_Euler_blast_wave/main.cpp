/**
 * @file main.cpp
 * @brief Implicit DG for 3D conservation laws
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-28
 */

#include <iostream>
#include <string>
#include <sstream>
#include "DGFEMSpace1D_GSL.h"

VEC<double> f0(const VEC<double>& u, double x, double t) {
  VEC<double> U(DIM);
  double rho(1), v(0), p(1e-9), pc(1e4);
  //if(x >= 0.495 && x <= 0.505) {
  if(x >= 0.5 && x <= 0.505) {
    U[0] = rho;
    U[1] = rho*v;
    U[2] = 0.5*rho*v*v + pc/(GAMMA-1.);
  }
  else {
    U[0] = rho;
    U[1] = rho*v;
    U[2] = 0.5*rho*v*v + p/(GAMMA-1.);
  }
  return U;
}

VEC<double> f(const VEC<double>& u) {
  VEC<double> F(DIM);
  double v = u[1]/u[0];
  double p = (GAMMA-1)*(u[2] - 0.5*u[0]*v*v);
  F[0] = u[1];
  F[1] = u[0]*v*v + p;
  F[2] = (u[2]+p)*v;
  return F;
}

VEC<double> source(const VEC<double>& u, double x, double t) {
  VEC<double> F(DIM);
  F[0] = 0.;
  F[1] = 0.;
  F[2] = 0.;
  return F;
}

VEC<VEC<double> > f_prime(const VEC<double>& u) {
  VEC<VEC<double> > a(DIM);
  for(u_int i = 0; i < DIM; ++i)
    a[i].resize(DIM, 0);
  double v = u[1]/u[0];
  double p = (GAMMA-1)*(u[2] - 0.5*u[0]*v*v);
  double H = (u[2]+p)/u[0];
  a[0][0] = 0; a[0][1] = 1; a[0][2] = 0;
  a[1][0] = 0.5*(GAMMA-3)*v*v; a[1][1] = (3-GAMMA)*v; a[1][2] = GAMMA-1;
  a[2][0] = v*(0.5*(GAMMA-1)*v*v-H); a[2][1] = H-(GAMMA-1)*v*v; a[2][2] = GAMMA*v;
  return a;
}

int main(int argc, char *argv[]) {
  if(argc != 5) {
    std::cout << "Usage: <Nx> <xl> <xr> <t_end> " << std::endl;
    abort();
  }

  clock_t t1, t2;
  u_int Nx = atoi(argv[1]);
  double xl = atof(argv[2]);
  double xr = atof(argv[3]);
  double t_end = atof(argv[4]);
  std::cout << "Set up problem..." << std::endl;
  DGFEMSpace1D Problem(Nx, xl, xr);
  std::cout << "Build quadrature info..." << std::endl;
  Problem.BuildQuad(K+1);
  std::cout << "Initialize..." << std::endl;
  Problem.init(f0);
  std::cout << "Start to solve..." << std::endl;
  t1 = clock();
  Problem.run_unsteady(f, f_prime, source, t_end);
  t2 = clock();
  std::stringstream s;
  s << "example8_Nx" << Nx << "_K" << K << "_PP" << PP_limiter << ".dat";
  std::string filename(s.str());
  std::ofstream out(filename.c_str());
  std::cout << "Print solution to " << filename << "..." << std::endl;
  Problem.print_solution(out);
  out.close();
  std::cout << "Time consumed: "
    //<< std::setiosflags(std::ios::scientific)
    << (t2-t1)/(double)CLOCKS_PER_SEC << std::endl;

  return 0;
}

