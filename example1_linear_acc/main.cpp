/**
 * @file main.cpp
 * @brief Implicit DG for 1D conservation laws
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include <iostream>
#include <string>
#include <sstream>
//#include "DGFEMSpace1D.h"
#include "DGFEMSpace1D_GSL.h"

VEC<double> f0(const VEC<double>& u, double x, double t) {
  VEC<double> U(DIM);
  U[0] = pow(sin(x), 2);
  return U;
}

VEC<double> f(const VEC<double>& u) {
  VEC<double> F(DIM);
  F[0] = u[0];
  return F;
}

VEC<double> source(const VEC<double>& u, double x, double t) {
  VEC<double> F(DIM);
  F[0] = pow(sin(x), 4);
  return F;
}

VEC<VEC<double> > f_prime(const VEC<double>& u) {
  VEC<VEC<double> > a(DIM);
  for(u_int i = 0; i < DIM; ++i)
    a[i].resize(DIM, 0);
  a[0][0] = 1;
  return a;
}

VEC<double> exact(const double x, const double t) {
  VEC<double> U(DIM);
  U[0] = (sin(4*x)-8*sin(2*x)+12*x)/32;
  return U;
}

int main(int argc, char *argv[]) {
  if(argc != 4) {
    std::cout << "Usage: <Nx> <xl> <xr> " << std::endl;
    abort();
  }

  clock_t t1, t2;
  u_int Nx = atoi(argv[1]);
  double xl = atof(argv[2]);
  double xr = atof(argv[3]);
  std::cout << "Set up problem..." << std::endl;
  DGFEMSpace1D Problem(Nx, xl, xr);
  std::cout << "Build quadrature info..." << std::endl;
  Problem.BuildQuad(K+1);
  std::cout << "Initialize..." << std::endl;
  Problem.init(f0);
  std::cout << "Start to solve..." << std::endl;
  t1 = clock();
  Problem.run_steady(f, f_prime, source);
  t2 = clock();
  std::stringstream s;
  s << "example1_Nx" << Nx << "_K" << K << "_PP" << PP_limiter << ".dat";
  std::string filename(s.str());
  std::ofstream out(filename.c_str());
  std::cout << "Print solution to " << filename << "..." << std::endl;
  Problem.print_solution_integral(out);
  out.close();
  std::cout << "Time consumed: "
    //<< std::setiosflags(std::ios::scientific)
    << (t2-t1)/(double)CLOCKS_PER_SEC << std::endl;

  //exact solution
  //std::ofstream out1("ex1_exact.dat");
	//out1.precision(8);
	//out1 << std::showpos;
  //out1.setf(std::ios::scientific);
  //for(u_int i = 0; i < 1e3; ++i) {
    //double center = (xr-xl)/1e3*(i+0.5)+xl;
    //out1 << center << " " << exact(center, 0) << "\n";
  //}
  //out1 << std::endl;
  //out1 << std::defaultfloat;
  //out1.close();

  //Eigen::setNbThreads(4);
  //std::cout << Eigen::nbThreads() << std::endl;
  return 0;
}

