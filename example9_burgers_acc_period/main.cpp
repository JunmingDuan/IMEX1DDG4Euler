/**
 * @file main.cpp
 * @brief Implicit DG for 1D conservation laws
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-12-11
 */

#include <iostream>
#include <string>
#include <sstream>
//#include "DGFEMSpace1D.h"
#include "DGFEMSpace1D_GSL.h"

VEC<double> f0(const VEC<double>& u, double x, double t) {
  VEC<double> U(DIM);
  U[0] = 1.0+sin(x);
  return U;
}

VEC<double> f(const VEC<double>& u) {
  VEC<double> F(DIM);
  F[0] = 0.5*pow(u[0],2);
  return F;
}

VEC<double> source(const VEC<double>& u, double x, double t) {
  VEC<double> F(DIM);
  F[0] = 0;
  return F;
}

VEC<VEC<double> > f_prime(const VEC<double>& u) {
  VEC<VEC<double> > a(DIM);
  for(u_int i = 0; i < DIM; ++i)
    a[i].resize(DIM, 0);
  a[0][0] = u[0];
  return a;
}

VEC<double> exact(const double x, const double t) {
  VEC<double> U(DIM);
  U[0] = sqrt(8-8*cos(x/4.));
  return U;
}

int main(int argc, char *argv[]) {
  if(argc != 5) {
    std::cout << "Usage: <Nx> <xl> <xr> " << std::endl;
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
  s << "example9_Nx" << Nx << "_K" << K << "_PP" << PP_limiter << ".dat";
  std::string filename(s.str());
  std::ofstream out(filename.c_str());
  std::cout << "Print solution to " << filename << "..." << std::endl;
  Problem.print_solution_integral(out);
  //Problem.print_solution_average(out);
  out.close();
  std::cout << "Time consumed: "
    << (t2-t1)/(double)CLOCKS_PER_SEC << std::endl;

  return 0;
}

