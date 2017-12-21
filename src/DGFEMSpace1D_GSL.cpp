/**
 * @file DGFEMSpace1D_GSL.cpp
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-19
 */

#include <stdio.h>
#include <stdlib.h>
#include "DGFEMSpace1D_GSL.h"

DGFEMSpace1D::DGFEMSpace1D(u_int Nx, double xl, double xr)
  : Nx(Nx), xl(xl), xr(xr) {
  mesh.resize(Nx+1);
  h = (xr - xl)/Nx;
  for(u_int i = 0; i < Nx+1; ++i) {
    mesh[i] = i*h;
  }
  sol.resize(Nx);
  sol1.resize(Nx);
  cell_average.resize(Nx);
  cell_val.resize(Nx);
  for(u_int i = 0; i < Nx; ++i) {
    sol[i].resize(K);
    sol1[i].resize(K);
    cell_average[i].resize(DIM);
    for(u_int k = 0; k < K; ++k) {
      sol[i][k].resize(DIM);
      sol1[i][k].resize(DIM);
    }
  }
  A = gsl_spmatrix_alloc(Nx*K*DIM, Nx*K*DIM);
  B = gsl_matrix_alloc(Nx*K*DIM, Nx*K*DIM);
  //p_full_solver = gsl_permutation_alloc(Nx*K*DIM);
  rhs = gsl_vector_alloc(Nx*K*DIM);
  vec_u1 = gsl_vector_alloc(Nx*K*DIM);
  vec_u2 = gsl_vector_alloc(Nx*K*DIM);
}

void DGFEMSpace1D::BuildQuad(u_int np) {
  std::vector<std::vector<double> > pw;
  pw = QUADINFO[0].LGL(np);
  TemQuad.set_np(np);
  TemQuad.set_points(pw[0]);
  TemQuad.set_weight(pw[1]);
  QUADINFO.resize(Nx);
  std::vector<double> gv(2), lv(2);
  lv[0] = -1, lv[1] = 1;
  for(u_int i = 0; i < Nx; ++i) {
    //init cell_val
    cell_val[i].resize(np);
    QUADINFO[i].set_np(np);
    gv[0] = mesh[i], gv[1] = mesh[i+1];
    QUADINFO[i].set_weight(pw[1]);
    std::vector<double> x(np);
    for(u_int j = 0; j < np; ++j) {
      local_to_global(pw[0][j], lv, gv, &x[j]);
      cell_val[i][j].resize(DIM);
    }
    QUADINFO[i].set_points(x);
    QUADINFO[i].set_jacobi( local_to_global_jacobian(lv, gv),
        global_to_local_jacobian(lv, gv) );
  }
}

void DGFEMSpace1D::Projection(u_int cell, func f0, double t, bU& u) {
  //set to zero
  for(u_int k = 0; k < K; ++k) {
    for(u_int d = 0; d < DIM; ++d) {
      u[k][d] = 0;
    }
  }
  std::vector<double> x = TemQuad.points();
  std::vector<double> p = QUADINFO[cell].points();
  std::vector<double> w = QUADINFO[cell].weight();
  VEC<double> U(DIM), u0(DIM);
  std::vector<double> V;
  double jab = QUADINFO[cell].l2g_jacobian();
  for(u_int k = 0; k < K; ++k) {
    double basis(0);
    for(u_int g = 0; g < x.size(); ++g) {
      U = f0(u0, p[g], t);
      V = Poly(x[g]);
      u[k] += jab*U*V[k]*w[g];
      basis += jab*V[k]*V[k]*w[g];
    }
    u[k] /= basis;
  }
}

VEC<double> DGFEMSpace1D::Composition(const SOL& sol, u_int cell, double x, double t) {
  VEC<double> u(DIM);
  std::vector<double> lv(2), gv(2), V;
  lv[0] = -1, lv[1] = 1;
  gv[0] = mesh[cell], gv[1] = mesh[cell+1];
  double lp;
  global_to_local(x, lv, gv, &lp);
  V = Poly(lp);
  for(u_int k = 0; k < K; ++k) {
    u += sol[cell][k]*V[k];
  }

  return u;
}

void DGFEMSpace1D::Pk2val(const SOL& sol, VEC<VEC<VEC<double>>>& val) {
  std::vector<double> p;
  for(u_int i = 0; i < Nx; ++i) {
    p = QUADINFO[i].points();
    for(u_int g = 0; g < p.size(); ++g) {
      val[i][g] = Composition(sol, i, p[g], 0);
    }
  }
}

void DGFEMSpace1D::init(func f0) {
  for(u_int i = 0; i < Nx; ++i) {
    Projection(i, f0, 0, sol[i]);
  }
  //print_solution(std::cout);
}

double DGFEMSpace1D::cal_dt(const SOL& sol, afunc g) {
  //return 10*h;
  //return 0.266*h;//ex4
  //return 0.5*h/cal_characteristic_speed(sol, g);//ex5
  //return 0.5*h/cal_characteristic_speed(sol, g);//ex6
  //return 2*h/cal_characteristic_speed(sol, g);//ex7
  //return 2*h/cal_characteristic_speed(sol, g);//ex8
  return h/cal_characteristic_speed(sol, g);//ex9
}

double DGFEMSpace1D::cal_characteristic_speed(const SOL& sol, afunc g) {
  double center(0), a(0);
  VEC<double> u;
  VEC<VEC<double>> tmp;
  for(u_int i = 0; i < Nx; ++i) {
    if(DIM == 1) {
      center = 0.5*(mesh[i]+mesh[i+1]);
      tmp = g(Composition(sol,i,center,0));
      if(fabs(tmp[0][0]) > a) a = fabs(tmp[0][0]);
    }
    else if(DIM == 3) {//Euler, return abs(u)+c
      center = 0.5*(mesh[i]+mesh[i+1]);
      u = Composition(sol,i,center,0);
      double v = u[1]/u[0];
      double p = (GAMMA-1)*(u[2] - 0.5*u[0]*v*v);
      double c = sqrt(GAMMA*p/u[0]);
      if(c+fabs(v) > a) a = c+fabs(v);
      //std::cout << "acoutic speed, c: " << c  << " , i: " << i << std::endl;
      //std::cout << u << std::endl;
      //std::cout << center << std::endl;
    }
  }
  return a;
}

int DGFEMSpace1D::forward_one_step(const SOL& sol, const F FLUX, afunc g, func source,
    double t, double dt, double* dtt, SOL& sol_new) {
  double alpha = cal_characteristic_speed(sol, g);
  Newton_iter(sol, FLUX, g, source, t, dt, alpha, sol_new);
  return 0;
}

VEC<double> DGFEMSpace1D::LxF(const F FLUX, const VEC<double>& a, const VEC<double>& b, const double alpha) {
  VEC<double> lf(a.size());
  lf = 0.5*(FLUX(a)+FLUX(b) - alpha*(b-a));
  return lf;
}

void DGFEMSpace1D::Newton_iter(const SOL& sol, const F FLUX, const afunc g, func source,
    const double t, const double dt, const double alpha, SOL& sol_new) {
  int Nt_ite(0);
  double Nt_err(1), Fval_norm(1);
  sol_new = sol;
  SOL2EVEC(sol_new, vec_u1);
  //while (Nt_err > Nt_tol && Fval_norm > Nt_Ftol && Nt_ite < MaxNt_ite) {
  while (Nt_err > Nt_tol && Nt_ite < MaxNt_ite) {
    form_jacobian_rhs(sol_new, sol, FLUX, g, source, t, dt, alpha);
    //std::cout << "=====sol^n,A,rhs,sol^{n+1}=====" << std::endl;
    //std::cout << "vec_u1:" << std::endl;
    //gsl_vector_fprintf(stdout, vec_u1, "%.6lf");
    //std::cout << "A:" << std::endl;
    //gsl_spmatrix_fprintf(stdout, A, "%.6lf");
    //std::cout << "rhs:" << std::endl;
    //gsl_vector_fprintf(stdout, rhs, "%.6lf");
    //std::cout << "rhs_norm:" << std::endl;
    //std::cout << gsl_blas_dnrm2(rhs) << std::endl;
    solve_leqn(A, rhs, vec_u2);
    Nt_err = gsl_blas_dnrm2(vec_u2);
    gsl_vector_add(vec_u1, vec_u2);
    //std::cout << "vec_u1:" << std::endl;
    //gsl_vector_fprintf(stdout, vec_u1, "%.6lf");
    EVEC2SOL(sol_new, vec_u1);
    //gsl_vector * tmp = gsl_vector_alloc(Nx*K*DIM);
    //*tmp = NLF(FLUX, sol_new, sol, source, alpha, t, dt);
    //Fval_norm = gsl_blas_dnrm2(tmp);
    std::cout << "Nt_ite: " << ++Nt_ite
      << ", Nt_err: " << Nt_err
      //<< ", Fval: " << Fval_norm
      //<< ", residual2: " << cal_err(sol_new, 2)[0]
      << std::endl;
  }
}

int kronecker(const int a, const int b) {
  if(a == b) return 1;
  else return 0;
}

EVEC DGFEMSpace1D::NLF(const F FLUX, const SOL& sol, const SOL& soln, func source,
    const double alpha, const double t, const double dt) {
  EVEC * fk = gsl_vector_alloc(Nx*K*DIM);
  int row;
  std::vector<double> x = TemQuad.points();
  std::vector<double> pnt, wei;
  std::vector<double> PolyVal, PGVal;
  std::vector<double> lv(2), gv(2);
  lv[0] = -1, lv[1] = 1;
  gsl_vector_set_zero(fk);
  bU tmp_u(K);
  for(u_int k = 0; k < K; ++k) {
    tmp_u[k].resize(DIM);
  }
  for(u_int i = 0; i < Nx; ++i) {
    pnt = QUADINFO[i].points();
    wei = QUADINFO[i].weight();
    Projection(i, source, t, tmp_u);//projection of the source term
    //std::cout << "i: " << i << std::endl;
    for(u_int k = 0; k < K; ++k) {
      VEC<double> fu(DIM), fu1(DIM);
      VEC<double> UL(DIM), UR(DIM);
      for(u_int d = 0; d < DIM; ++d) {
        row = i*(K*DIM) + k*DIM + d;//fixed row
        //time derivative
        gsl_vector_set(fk, row, (sol[i][k][d]-soln[i][k][d])/(2*k+1));
        //std::cout << "time derivative:\n" << "row: " << row << ", fk: " << *(fk->data+row) << "\n";
        //element integral
        for(u_int g = 0; g < pnt.size(); ++g) {
          PGVal = PolyG(x[g]);
          fu = FLUX(Composition(sol,i,pnt[g],t));
          *(fk->data+row) += - dt/2 * fu[d] * (PGVal[k]*QUADINFO[i].g2l_jacobian()) * wei[g];
        }
        //std::cout << "element integral:\n" << "row: " << row << ", fk: " << *(fk->data+row) << "\n";
        //
        gv[0] = mesh[i], gv[1] = mesh[i+1];
        VEC<double> flux;
        int ghost_i;
        double val, ghost_local_pnt;

        //flux on the right boundary of the cell
        PolyVal = Poly(lv[1]);
        UL = Composition(sol,i,gv[1],t);
        if(i < Nx-1) {
          UR = Composition(sol,i+1,gv[1],t);
        }
        else if(i == Nx-1) {
          UR = RightBoundary(sol, UL, &ghost_i, &ghost_local_pnt, t);
        }
        flux = LxF(FLUX, UL, UR, alpha);
        val = flux[d] * dt/(gv[1]-gv[0]) * PolyVal[k];
        //std::cout << "outer:\n" << "row: " << row << ", fk: " << *(fk->data+row) << "\n";
        *(fk->data+row) += val;

        //flux on the left boundary of the cell
        PolyVal = Poly(lv[0]);
        UR = Composition(sol,i,gv[0],t);
        if(i > 0) {
          UL = Composition(sol,i-1,gv[0],t);
        }
        else if(i == 0) {
          UL = LeftBoundary(sol, UR, &ghost_i, &ghost_local_pnt, t);
        }
        flux = LxF(FLUX, UL, UR, alpha);
        val = flux[d] * dt/(gv[1]-gv[0]) * PolyVal[k];
        *(fk->data+row) -= val;
        //std::cout << "inner:\n" << "row: " << row << ", fk: " << *(fk->data+row) << "\n";

        //source term
        *(fk->data+row) -= tmp_u[k][d] * dt/(2*k+1);
        //std::cout << "source:\n" << "row: " << row << ", fk: " << *(fk->data+row) << "\n";
      }
    }
  }
  return *fk;
}

void my_gsl_spmat_add(gsl_spmatrix*A, const int row, const int col, const double val) {
  double * ptr = gsl_spmatrix_ptr(A, row, col);
  if(ptr == NULL) gsl_spmatrix_set(A, row, col, val);
  else *ptr += val;
}

/**
 * @brief form_jacobian_rhs
 *        jacobian matrix is a (Nx*K) dimensional square matrix.
 *        rhs is a (Nx*K) dimensional vector, we order the vector
 *        in space firstly, then in polynomial space and last in physical space.
 *        rhs[i*(K*DIM)+k*DIM+d], is the d-th physical variable of the k-th
 *        polynomial in the i-the cell.
 *
 * @param sol current sol
 *        soln sol at t^n
 */
void DGFEMSpace1D::form_jacobian_rhs(const SOL& sol, const SOL& soln, const F FLUX, afunc fp, func source,
    const double t, const double dt, const double alpha) {
  int row, col;
  double val, val0;
  std::vector<double> pnt, wei;
  std::vector<double> x = TemQuad.points();
  std::vector<double> PolyVal, PGVal, LocPolyVal;
  std::vector<double> lv(2), gv(2);
  lv[0] = -1, lv[1] = 1;
  gsl_spmatrix_set_zero(A);
  for(u_int i = 0; i < Nx; ++i) {//i-th cell
    //std::cout << "i: " << i << std::endl;
    pnt = QUADINFO[i].points();
    wei = QUADINFO[i].weight();
    for(u_int k = 0; k < K; ++k) {//k-th phi
      //std::cout << "k: " << k << std::endl;
      for(u_int d = 0; d < DIM; ++d) {//d-th equation
        //std::cout << "d: " << d << std::endl;
        row = i*(K*DIM) + k*DIM + d;//fixed row
        //time derivative: u_{i,d}^(k)
        col = row;
        val = 1./(2*k+1);
        my_gsl_spmat_add(A, row, col, val);
        //std::cout << "time, derivative:\n" << "col: " << col << ", val: " << val << "\n";

        VEC<VEC<double> > AF;
        VEC<double> UL(DIM), UR(DIM);
        //element integral: u^(q)_{i,d}, q=0..K-1
        for(u_int q = 0; q < K; ++q) {//derivative of u^(q)_{i,m}
          //std::cout << "q: " << q << std::endl;
          for(u_int m = 0; m < DIM; ++m) {
            //std::cout << "m: " << m << std::endl;
            val = 0;
            col = i*(K*DIM) + q*DIM + m;
            for(u_int g = 0; g < pnt.size(); ++g) {//numerical integral
              //std::cout << "g: " << g << std::endl;
              VEC<double> U0 = Composition(sol,i, pnt[g], 0);
              AF = fp(U0);
              PolyVal = Poly(x[g]);
              PGVal = PolyG(x[g]);
              //only the m-th component of U, i.e.,
              //U_{i,m}=\sum u^(q)_{i,m}*poly^(q) has component u^(q)_{i,m}
              double pd = AF[d][m] * PolyVal[q];
              //scale in the prime of the polynomial
              val += pd * (PGVal[k]*QUADINFO[i].g2l_jacobian()) * wei[g];
            }
            val *= -dt/2;
            my_gsl_spmat_add(A, row, col, val);
            //std::cout << "element integral:\n" << "col: " << col << ", val: " << val << "\n";
          }
        }
        //flux on the outer boundary
        gv[0] = mesh[i], gv[1] = mesh[i+1];
        PolyVal = Poly(lv[1]);
        int ghost_i;
        double ghost_local_pnt;
        for(u_int q = 0; q < K; ++q) {//derivative of u^(q)_{i,m} and u^(q)_{i+1,m}
          for(u_int m = 0; m < DIM; ++m) {
            //i-th cell
            LocPolyVal = Poly(lv[1]);
            col = i*(K*DIM) + q*DIM + m;
            VEC<double> U0 = Composition(sol,i, gv[1], 0);
            AF = fp(U0);
            double pd = 0.5 * (AF[d][m]+alpha*kronecker(d,m)) * LocPolyVal[q];
            val = pd * dt/(gv[1]-gv[0]) * PolyVal[k];
            my_gsl_spmat_add(A, row, col, val);
            //std::cout << "outer, ith:\n" << "col: " << col << ", val: " << val << "\n";

            //(i+1)-th cell
            if(i < Nx-1) {
              LocPolyVal = Poly(lv[0]);
              col = (i+1)*(K*DIM) + q*DIM + m;
              UR = Composition(sol,i+1, gv[1], 0);
            }
            else if(i == Nx-1) {
              UR = RightBoundary(sol, U0, &ghost_i, &ghost_local_pnt, t);
              LocPolyVal = Poly(ghost_local_pnt);
              col = (ghost_i)*(K*DIM) + q*DIM + m;//modified here
            }
            AF = fp(UR);
            pd = 0.5 * (AF[d][m]-alpha*kronecker(d,m)) * LocPolyVal[q];
            val = pd * dt/(gv[1]-gv[0]) * PolyVal[k];
            my_gsl_spmat_add(A, row, col, val);
            //std::cout << "outer (i+1)th:\n" << "col: " << col << ", val: " << val << "\n";
          }
        }

        //flux on the inner boundary
        PolyVal = Poly(lv[0]);
        for(u_int q = 0; q < K; ++q) {//derivative of u^(q)_{i,m} and u^(q)_{i-1,m}
          for(u_int m = 0; m < DIM; ++m) {
            double pd;
            //i-th cell
            LocPolyVal = Poly(lv[0]);
            col = i*(K*DIM) + q*DIM + m;
            VEC<double> U0 = Composition(sol,i, gv[0], 0);
            AF = fp(U0);
            pd = 0.5 * (AF[d][m]-alpha*kronecker(d,m)) * LocPolyVal[q];
            val = - pd * dt/(gv[1]-gv[0]) * PolyVal[k];
            my_gsl_spmat_add(A, row, col, val);
            //std::cout << "inner (i)th:\n" << "col: " << col << ", val: " << val << "\n";

            //(i-1)-th cell
            if(i > 0) {
              LocPolyVal = Poly(lv[1]);
              col = (i-1)*(K*DIM) + q*DIM + m;
              UL = Composition(sol,i-1, gv[0], 0);
            }
            else if(i == 0) {
              UL = LeftBoundary(sol, U0, &ghost_i, &ghost_local_pnt, t);
              LocPolyVal = Poly(ghost_local_pnt);
              col = (ghost_i)*(K*DIM) + q*DIM + m;
            }
            AF= fp(UL);
            pd = 0.5 * (AF[d][m]+alpha*kronecker(d,m)) * LocPolyVal[q];
            val = - pd * dt/(gv[1]-gv[0]) * PolyVal[k];
            if(i == 0 && BDL == 0) val = 0;//ghost = 0
            my_gsl_spmat_add(A, row, col, val);
            //std::cout << "inner (i-1)th:\n" << "col: " << col << ", val: " << val << "\n";
          }
        }
        //gsl_spmatrix_fprintf(stdout, A, "%.6e");
        //abort();

      }
    }
  }
  //RHS
  *rhs = NLF(FLUX, sol, soln, source, alpha, t, dt);
  gsl_vector_scale(rhs, -1);

  //std::cout << "A:\n" << std::endl;
  //gsl_spmatrix_fprintf(stdout, A, "%.6e");
  //std::cout << "B:\n" << std::endl;
  //for(size_t i = 0; i < Nx*K*DIM; ++i) {
    //for(size_t j = 0; j < Nx*K*DIM; ++j) {
      //printf("%+.6e ", gsl_matrix_get(B,i,j));
    //}
    //printf("\n");
  //}
  //gsl_matrix_fprintf(stdout, B, "%.6e");
  //std::cout.precision(16);
  //std::cout << std::showpos;
  //std::cout.setf(std::ios::scientific);
  //for(u_int i = 0; i < Nx*K*DIM; ++i) {
    //for(u_int j = 0; j < Nx*K*DIM; ++j) {
      //std::cout << gsl_spmatrix_get(A,i,j) << " ";
    //}
    //std::cout << "\n";
  //}
  //std::cout << std::endl;
  //gsl_vector_fprintf(stdout, rhs, "%+.6e");
  //std::cout << std::defaultfloat;
  //abort();
}

void DGFEMSpace1D::solve_leqn(MAT *A, const EVEC *rhs, EVEC *u) {
  //std::cout << "======solve_leqn by GSL GMRES======" << std::endl;
  const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
  gsl_splinalg_itersolve *work =
    gsl_splinalg_itersolve_alloc(T, Nx*K*DIM, 0);
  int iter = 0;
  double residual;
  int status;
  /* initial guess u = 0 */
  gsl_vector_set_zero(u);
  /* solve the system A u = f */
  do
  {
    status = gsl_splinalg_itersolve_iterate(A, rhs, tol, u, work);
    /* print out residual norm ||A*u - f|| */
    residual = gsl_splinalg_itersolve_normr(work);
    //fprintf(stdout, "iter %zu residual = %.12e\n", iter++, residual);

    //if (status == GSL_SUCCESS)
      //fprintf(stdout, "Converged, iter %d, residual %.6e\n", ++iter, residual);
  }
  while (status == GSL_CONTINUE && iter < 1e1);
  //std::cout << "===================================" << std::endl;
  //std::cout << "==solve_leqn by GSL full_matrix_LU==" << std::endl;
  //int s;
  //gsl_spmatrix_sp2d(B, A);
  //gsl_permutation * p = gsl_permutation_alloc(Nx*K*DIM);
  //gsl_linalg_LU_decomp(B, p, &s);
  //gsl_linalg_LU_solve(B, p, rhs, u);
  //gsl_permutation_free(p);
  //std::cout << "====================================" << std::endl;
}

VEC<double> exact0(const VEC<double>& u, double x, double t) {
  VEC<double> U(DIM);
  //U[0] = (sin(4*x)-8*sin(2*x)+12*x)/32;
  //U[0] = 0;
  double s = cos(x/4.);
  U[0] = sqrt(8.*s*(s*s-3)/3+16./3);
  return U;
}

void DGFEMSpace1D::run_steady(F FLUX, afunc g, func source) {
  int ite(0), is_pp(0);
  double t(0), dt(0), dtt(0);
  VEC<double> err(DIM,1), tol(DIM,TOL);
  while ( err[0] > tol[0] ) {//only for 1D case
    dt = cal_dt(sol, g);
    if(PP_limiter == 1) {
      do {
        forward_one_step(sol, FLUX, g, source, t, dt, &dtt, sol1);
        is_pp = judge_positivity(sol1);
        if(is_pp == 0) dt *= 0.5;
      } while(is_pp == 0);
      Pk2val(sol1, cell_val);
      scaling_limiter::run(cell_average, cell_val, sol1);
    }
    else if(PP_limiter == 0) {
      forward_one_step(sol, FLUX, g, source, t, dt, &dtt, sol1);
    }
    err = cal_norm(sol, sol1, 2);
    sol = sol1;
    ite++;
    t += dt;
    std::cout << "ite: " << ite << ", dt: " << dt << ", t: " << t << ", err: ";
    std::cout << err << std::endl;
  }
}

void DGFEMSpace1D::run_unsteady(F FLUX, afunc g, func source, double t_end) {
  int ite(0), is_pp(0), circle(0);
  double t(0), dt(0), dtt(0);
  VEC<double> err(DIM,1), tol(DIM,TOL);
  while ( t < t_end ) {
    dt = cal_dt(sol, g);
    if(t + dt > t_end) dt = t_end - t;
    if(PP_limiter == 1) {
      circle = 0;
      do {
        forward_one_step(sol, FLUX, g, source, t, dt, &dtt, sol1);
        is_pp = judge_positivity(sol1);
        if(is_pp == 0) {
          dt *= 0.5;
        }
        circle++;
      } while(is_pp == 0);
      std::cout << "circle: " << circle << std::endl;
      Pk2val(sol1, cell_val);
      scaling_limiter::run(cell_average, cell_val, sol1);
    }
    else if(PP_limiter == 0) {
      forward_one_step(sol, FLUX, g, source, t, dt, &dtt, sol1);
    }
    err = cal_norm(sol, sol1, 2);
    sol = sol1;
    ite++;
    t += dt;
    std::cout << "ite: " << ite << ", dt: " << dt << ", t: " << t << std::endl;
  }
}


int DGFEMSpace1D::judge_positivity(const SOL& sol) {
  for(u_int i = 0; i < Nx; ++i) {
    //cell_average is the first moment
    cell_average[i] = sol[i][0];
    if(cell_average[i][0] < EPS) {
      std::cout << "av[i] < 0: " << i << " " << cell_average[i][0] << std::endl;
      return 0;
    }
    if(DIM == 3) {
      double e = cell_average[i][2]-0.5*pow(cell_average[i][1],2)/cell_average[i][0];
      if(e < EPS) {
        std::cout << "e[i] < 0: " << i << " " << e << std::endl;
        return 0;
      }
    }
  }
  return 1;
}

void DGFEMSpace1D::SOL2EVEC(const SOL& sol, EVEC *vec_u) {
  for(u_int i = 0; i < Nx; ++i) {
    for(u_int k = 0; k < K; ++k) {
      for(u_int d = 0; d < DIM; ++d) {
        gsl_vector_set(vec_u, i*(K*DIM)+k*DIM+d, sol[i][k][d]);
      }
    }
  }
}

void DGFEMSpace1D::EVEC2SOL(SOL& sol, const EVEC *vec_u) {
  for(u_int i = 0; i < Nx; ++i) {
    for(u_int k = 0; k < K; ++k) {
      for(u_int d = 0; d < DIM; ++d) {
        sol[i][k][d] = gsl_vector_get(vec_u, i*(K*DIM)+k*DIM+d);
      }
    }
  }
}

VEC<double> DGFEMSpace1D::cal_norm(const SOL& s1, const SOL& s2, int n) {
  VEC<double> norm(DIM), tmp1, tmp2;
  if(n == 2) {
    for(u_int i = 0; i < Nx; ++i) {
      std::vector<double> p = QUADINFO[i].points();
      std::vector<double> w = QUADINFO[i].weight();
      for(u_int g = 0; g < p.size(); ++g) {
        tmp1 = Composition(s1,i,p[g],0);
        tmp2 = Composition(s2,i,p[g],0);
        for(u_int d = 0; d < DIM; ++d) {
          norm[d] += pow(tmp1[d]-tmp2[d], 2) * w[g];
        }
      }
    }
    for(u_int d = 0; d < DIM; ++d) {
      norm[d] = sqrt(norm[d]/Nx);
    }
    return norm;
  }
  else {
    std::cout << "Wrong norm choice!" << std::endl;
    return norm;
  }
}

VEC<double> exact1(const double x, const double t) {
  VEC<double> U(DIM);
  double s = cos(x/4.);
  //U[0] = sqrt(8-8*cos(x/4.));
  U[0] = sqrt(8.*s*(s*s-3)/3+16./3);
  return U;
}

VEC<double> DGFEMSpace1D::cal_err(const SOL& s1, int n) {
  VEC<double> norm(DIM,0), tmp1, tmp2;
  double center;
  if(n == 2) {
    for(u_int i = 0; i < Nx; ++i) {
      center = 0.5*(mesh[i]+mesh[i+1]);
      tmp1 = Composition(s1,i,center,0);
      tmp2 = exact1(center,0);
      for(u_int d = 0; d < DIM; ++d) {
        norm[d] += pow(tmp1[d]-tmp2[d], 2);
      }
    }
    for(u_int d = 0; d < DIM; ++d) {
      norm[d] = sqrt(norm[d]/Nx);
    }
    return norm;
  }
  else if(n == 1) {
    for(u_int i = 0; i < Nx; ++i) {
      center = 0.5*(mesh[i]+mesh[i+1]);
      tmp1 = Composition(s1,i,center,0);
      tmp2 = exact1(center,0);
      for(u_int d = 0; d < DIM; ++d) {
        norm[d] += fabs(tmp1[d]-tmp2[d]);
      }
    }
    for(u_int d = 0; d < DIM; ++d) {
      norm[d] = (norm[d]/Nx);
    }
    return norm;
  }
  else if(n == 0) {
    for(u_int i = 0; i < Nx; ++i) {
      center = 0.5*(mesh[i]+mesh[i+1]);
      tmp1 = Composition(s1,i,center,0);
      tmp2 = exact1(center,0);
      for(u_int d = 0; d < DIM; ++d) {
        if( fabs(tmp1[d]-tmp2[d]) > norm[d] ) norm[d] = fabs(tmp1[d]-tmp2[d]);
      }
    }
    return norm;

  }
  else {
    std::cout << "Wrong norm choice!" << std::endl;
    return norm;
  }
}

//VEC<double> DGFEMSpace1D::cal_err(const SOL& s1, int n) {
  //VEC<double> norm(DIM,0), tmp1, tmp2;
  //if(n == 2) {
    //for(u_int i = 0; i < Nx; ++i) {
      //std::vector<double> p = QUADINFO[i].points();
      //std::vector<double> w = QUADINFO[i].weight();
      //for(u_int g = 0; g < p.size(); ++g) {
        //tmp1 = Composition(s1,i,p[g],0);
        //tmp2 = exact1(p[g],0);
        //for(u_int d = 0; d < DIM; ++d) {
          //norm[d] += pow(tmp1[d]-tmp2[d], 2) * w[g];
        //}
      //}
    //}
    //for(u_int d = 0; d < DIM; ++d) {
      //norm[d] = sqrt(norm[d]/Nx/2);
    //}
    //return norm;
  //}
  //else if(n == 1) {
    //for(u_int i = 0; i < Nx; ++i) {
      //std::vector<double> p = QUADINFO[i].points();
      //std::vector<double> w = QUADINFO[i].weight();
      //for(u_int g = 0; g < p.size(); ++g) {
        //tmp1 = Composition(s1,i,p[g],0);
        //tmp2 = exact1(p[g],0);
        //for(u_int d = 0; d < DIM; ++d) {
          //norm[d] += fabs(tmp1[d]-tmp2[d]) * w[g];
        //}
      //}
    //}
    //for(u_int d = 0; d < DIM; ++d) {
      //norm[d] = (norm[d]/Nx/2);
    //}
    //return norm;
  //}
  //else if(n == 0) {
    //for(u_int i = 0; i < Nx; ++i) {
      //std::vector<double> p = QUADINFO[i].points();
      //std::vector<double> w = QUADINFO[i].weight();
      //for(u_int g = 0; g < p.size(); ++g) {
        //tmp1 = Composition(s1,i,p[g],0);
        //tmp2 = exact1(p[g],0);
        //for(u_int d = 0; d < DIM; ++d) {
          //if( fabs(tmp1[d]-tmp2[d]) > norm[d] ) norm[d] = fabs(tmp1[d]-tmp2[d]);
        //}
      //}
    //}
    //return norm;

  //}
  //else {
    //std::cout << "Wrong norm choice!" << std::endl;
    //return norm;
  //}
//}

//void DGFEMSpace1D::print_solution(std::ostream& os) {
	//os.precision(16);
	//os << std::showpos;
  //os.setf(std::ios::scientific);
  //double center;
  //for(u_int i = 0; i < Nx; ++i) {
    //center = 0.5*(mesh[i]+mesh[i+1]);
    //os << center << " "  << Composition(sol,i,center,0) << "\n";
  //}
  //os << std::endl;
  //os << std::defaultfloat;
//}

void DGFEMSpace1D::print_solution_integral(std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  std::vector<double> x = TemQuad.points();
  for(u_int i = 0; i < Nx; ++i) {
    std::vector<double> w = QUADINFO[i].weight();
    std::vector<double> p = QUADINFO[i].points();
    for(u_int g = 0; g < x.size(); ++g) {
      os << p[g] << " "  << w[g] << " " << Composition(sol,i,p[g],0) << "\n";
    }
    os << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

void DGFEMSpace1D::print_solution_average(std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  for(u_int i = 0; i < Nx; ++i) {
      os << 0.5*(mesh[i]+mesh[i+1]) << " " << sol[i][0] << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

MAT DGFEMSpace1D::get_A() {
  return *A;
}

