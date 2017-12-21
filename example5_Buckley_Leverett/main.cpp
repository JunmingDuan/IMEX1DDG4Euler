/**
 * @file main.cpp
 * @brief Implicit DG for 1D conservation laws
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-27
 */

#include <iostream>
#include <string>
#include <sstream>
#include "DGFEMSpace1D_GSL.h"

VEC<double> f0(const VEC<double>& u, double x, double t) {
  VEC<double> U(DIM);
  if(x >= 3. && x <= 4.) U[0] = 0.9;
  else U[0] = 1e-3;
  return U;
}

VEC<double> f(const VEC<double>& u) {
  VEC<double> F(DIM);
  double s = u[0];
  F[0] = 4.*s*s/(4*s*s+pow(1-s,2));
  return F;
}

VEC<double> source(const VEC<double>& u, double x, double t) {
  VEC<double> F(DIM);
  F[0] = 0.;
  return F;
}

VEC<VEC<double> > f_prime(const VEC<double>& u) {
  VEC<VEC<double> > a(DIM);
  for(u_int i = 0; i < DIM; ++i)
    a[i].resize(DIM, 0);
  double s = u[0];
  a[0][0] = 8*s/(4*s*s+pow(1-s,2)) - (4*s*s*(10*s-2))/pow(4*s*s+pow(1-s,2),2);
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
  double Nt_tol(1e-13), Nt_Ftol(1e-13), TOL(1e-6);
  std::cout << "Set up problem..." << std::endl;
  DGFEMSpace1D Problem(Nx, xl, xr, Nt_tol, Nt_Ftol, TOL);
  std::cout << "Build quadrature info..." << std::endl;
  Problem.BuildQuad(K+1);
  std::cout << "Initialize..." << std::endl;
  Problem.init(f0);
  std::cout << "Start to solve..." << std::endl;
  t1 = clock();
  Problem.run_unsteady(f, f_prime, source, t_end);
  t2 = clock();
  std::stringstream s;
  s << "example5_Nx" << Nx << "_K" << K << ".dat";
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

