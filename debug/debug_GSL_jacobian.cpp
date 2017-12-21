//#include <cstdlib>
#include "DGFEMSpace1D_GSL.h"

VEC<double> f0(const VEC<double>& u, double x, double t) {
  VEC<double> U(DIM);
  double rhol(1), rhor(1), ul(0), ur(0), pl(1e3), pr(1e-2);
  if(x < 0.5) {
    U[0] = rhol;
    U[1] = rhol*ul;
    U[2] = 0.5*rhol*ul*ul + pl/(GAMMA-1.);
  }
  else {
    U[0] = rhor;
    U[1] = rhor*ur;
    U[2] = 0.5*rhor*ur*ur + pr/(GAMMA-1.);
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

//debug jacobian by differential F
//verifying the jacobian matrix
int main(int argc, char* argv[]) {
  if(argc < 4) {
    std::cout << "Usage: <Nx> <xl> <xr> " << std::endl;
    abort();
  }
  u_int Nx = atoi(argv[1]);
  double xl = atof(argv[2]);
  double xr = atof(argv[3]);
  double Nt_tol(1e-14), Nt_Ftol(1e-14), TOL(1e-15);

  std::cout << "Set up problem..." << std::endl;
  DGFEMSpace1D Problem(Nx, xl, xr, Nt_tol, Nt_Ftol, TOL);
  std::cout << "Build quadrature info..." << std::endl;
  Problem.BuildQuad(K+1);
  std::cout << "Initialize..." << std::endl;
  Problem.init(f0);

  std::cout << "########## DEBUG jacobian #########" << std::endl;
  SOL x1, x2, xn;
  double eps(1e-10), a(0), dt(0.1);
  int i1, k1, d1, i2, k2, d2;
  for(u_int loop = 0; loop < 1e3; ++loop) {
    std::cin.get();
    //i1 = 2, k1 = 0, d1 = 2;
    //i2 = 1, k2 = 0, d2 = 0;
    srand((unsigned)time(NULL));
    a = rand()/(double)RAND_MAX;
    //a = 1;
    i1 = rand()%Nx, k1 = rand()%K, d1 = rand()%DIM;
    i2 = rand()%Nx, k2 = rand()%K, d2 = rand()%DIM;
    int row = i1*(K*DIM)+k1*DIM+d1;
    int col = i2*(K*DIM)+k2*DIM+d2;
    //initialize xn, x1, x2
    x1.resize(Nx); x2.resize(Nx); xn.resize(Nx);
    for(u_int i = 0; i < Nx; ++i) {
      x1[i].resize(K); x2[i].resize(K); xn[i].resize(K);
      for(u_int k = 0; k < K; ++k) {
        x1[i][k].resize(DIM); x2[i][k].resize(DIM); xn[i][k].resize(DIM);
        for(u_int d = 0; d < DIM; ++d) {
          xn[i][k][d] = rand()/(double)RAND_MAX;
          x2[i][k][d] = rand()/(double)RAND_MAX;
          x1[i][k][d] = x2[i][k][d];
          //xn[i][k][d] = 1;
          //x2[i][k][d] = 1;
          //x1[i][k][d] = 1;
        }
      }
    }
    x1[i2][k2][d2] = x2[i2][k2][d2]+eps;

    //std::cout << "r1:\n" << std::endl;
    EVEC r1 = Problem.NLF(f, x1, xn, source, a, 0, dt);
    //std::cout << "r2:\n" << std::endl;
    EVEC r2 = Problem.NLF(f, x2, xn, source, a, 0, dt);
    //std::cout << "jacobian:\n" << std::endl;
    Problem.form_jacobian_rhs(x2, xn, f, f_prime, source, 0, dt, a);
    MAT B = Problem.get_A();

    double tmp = (r1.data[row]-r2.data[row])/eps - gsl_spmatrix_get(&B, row, col);
    std::cout << "DIFF-A_ij: " << tmp << std::endl;
    //std::cout << "i1,k1,d1: " << i1 << " " << k1 << " " << d1 << std::endl;
    //std::cout << "i2,k2,d2: " << i2 << " " << k2 << " " << d2 << std::endl;
    std::cout << "row,col: " << row << " " << col << std::endl;
    std::cout << "DIFF[row]: " << (r1.data[row]-r2.data[row])/eps << " ,A_ij: " << gsl_spmatrix_get(&B, row, col) << std::endl;
    if(fabs(tmp) > 1e-4) {
      std::cout << "row,col: " << row << " " << col << std::endl;
      std::cout << "i1,k1,d1: " << i1 << " " << k1 << " " << d1 << std::endl;
      std::cout << "i2,k2,d2: " << i2 << " " << k2 << " " << d2 << std::endl;
      //std::cout << "NLF(x1, xn): " << Problem.NLF(f, x1, xn, a, 0, dt)[row]
        //<< " NLF(x2, xn): " << Problem.NLF(f, x2, xn, a, 0, dt)[row] << std::endl;
      //std::cout << "DIFF[row]: " << DIFF[row] << " ,A_ij: " << Problem.get_A().coeffRef(row, col) << std::endl;
      std::cout << "xn: " << xn << "\n";
      std::cout << "x1: " << x1 << "\n";
      std::cout << "x2: " << x2 << std::endl;
      abort();
    }
  }

  std::cout << "########## DEBUG jacobian ##########" << std::endl;

  return 0;
}

