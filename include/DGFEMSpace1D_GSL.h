/**
 * @file DGFEMSpace1D_GSL.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-19
 */

#ifndef DGFEMSPACE1D_H
#define DGFEMSPACE1D_H
#include "gsl/gsl_math.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_spmatrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_splinalg.h"
#include "BasFun.h"
#include "Quadrature.h"
#include "interval_crd_trs.h"
#include "scaling_limiter.h"

typedef int BM;
typedef std::vector<std::vector<double> > QUAD;
typedef VEC<double> (*func)(const VEC<double>&,double,double);
typedef VEC<VEC<double> > (*afunc)(const VEC<double>&);
typedef VEC<double> (*F)(const VEC<double>&);
typedef gsl_spmatrix MAT;
typedef gsl_vector EVEC;

class DGFEMSpace1D {
  private:
    u_int Nx;
    double xl, xr;
    double h;
    VEC<double> mesh;
    TemplateQuadrature TemQuad;
    std::vector<Quadrature> QUADINFO;
    SOL sol, sol1;//dimension: Nx*K*DIM
    VEC<VEC<double>> cell_average;//dimension: Nx*DIM
    VEC<VEC<VEC<double>>> cell_val;//dimension: Nx*G*DIM
    BM bml, bmr;
    MAT *A;
    gsl_matrix *B;
    EVEC *rhs, *vec_u1, *vec_u2;
    //gsl_splinalg_itersolve_type * solver;
    //gsl_permutation * p_full_solver;
    //gsl_linalg_LU_solve * full_solver;

  public:
    DGFEMSpace1D(u_int Nx, double xl, double xr);
    void BuildQuad(u_int np);
    void Projection(u_int cell, func f0, double t, bU&);
    VEC<double> Composition(const SOL&, u_int cell, double x, double t);
    void Pk2val(const SOL&, VEC<VEC<VEC<double>>>&);
    void init(func f0);
    double cal_dt(const SOL&, afunc);
    double cal_characteristic_speed(const SOL&, afunc);
    /**
     * @brief forward_one_step
     *
     * @param double dt
     *
     * @return 0, donnot change dt; 1, change dt to dtt
     */
    int forward_one_step(const SOL&, const F, afunc, func,
        double t, double dt, double* dtt, SOL&);
    /**
     * @brief Newton_iter
     *
     * @param sol solution at t^n
     * @param dt
     */
    VEC<double> LxF(const F, const VEC<double>&, const VEC<double>&, const double);
    void Newton_iter(const SOL& sol, const F, const afunc g, func,
        const double, const double, const double, SOL&);
    /**
     * @brief NLF nonlinear function
     *
     * @param sol
     *
     * @return RHS, i.e., F(sol)
     */
    VEC<double> LeftBoundary(const SOL& sol, const VEC<double>& U0, int*, double*, const double t);
    VEC<double> RightBoundary(const SOL& sol, const VEC<double>& U0, int*, double*, const double t);
    EVEC NLF(const F, const SOL& sol, const SOL& soln, func,
        const double alpha, const double t, const double dt);
    void form_jacobian_rhs(const SOL& sol, const SOL& soln, const F, afunc, func,
        const double, const double, const double);
    void solve_leqn(MAT*, const EVEC*, EVEC*);
    void run_steady(F, afunc, func);
    void run_unsteady(F, afunc, func, double t_end);
    int judge_positivity(const SOL&);
    void SOL2EVEC(const SOL&, EVEC*);
    void EVEC2SOL(SOL&, const EVEC*);
    VEC<double> cal_norm(const SOL&, const SOL&, int);
    VEC<double> cal_err(const SOL& s1, int n);
    void print_solution_integral(std::ostream&);
    void print_solution_average(std::ostream&);
    MAT get_A();
};

#endif //DGFEMSPACE1D_H

