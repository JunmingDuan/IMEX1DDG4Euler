/**
 * @file DGFEMSpace1D.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#ifndef DGFEMSPACE1D_H
#define DGFEMSPACE1D_H
#include "vvector.h"
#include "BasFun.h"
#include "Quadrature.h"
#include "interval_crd_trs.h"
#include "Eigen/Sparse"
#include "Eigen/SparseLU"
#include "Eigen/IterativeLinearSolvers"
#include "unsupported/Eigen/IterativeSolvers"

#define VEC vvector
typedef VEC<VEC<double> > bU;
typedef VEC<bU> SOL;
typedef int BM;
typedef std::vector<std::vector<double> > QUAD;
typedef VEC<double> (*func)(const VEC<double>&,double,double);
typedef VEC<VEC<double> > (*afunc)(const VEC<double>&);
typedef VEC<double> (*F)(const VEC<double>&);
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> MAT;
typedef Eigen::VectorXd EVEC;

/**
 * @brief dimension of the equation, 1 for scalar equation and 3 for Euler equations
 */
const u_int DIM = 1;

class DGFEMSpace1D {
  private:
    u_int Nx;
    double xl, xr;
    double h;
    double Nt_tol, Nt_Ftol, TOL;
    VEC<double> mesh;
    TemplateQuadrature TemQuad;
    std::vector<Quadrature> QUADINFO;
    SOL sol, sol1;
    BM bml, bmr;
    MAT A;
    EVEC rhs, vec_u1, vec_u2;
    Eigen::SparseLU<MAT> solver;
    //Eigen::SparseLU<MAT, Eigen::NaturalOrdering<int>> solver;
    //Eigen::BiCGSTAB<MAT> solver;
    //Eigen::BiCGSTAB<MAT, Eigen::IncompleteLUT<double>> solver;
    //Eigen::BiCGSTAB<MAT, Eigen::IdentityPreconditioner> solver;
    //Eigen::BiCGSTAB<MAT, Eigen::DiagonalPreconditioner<double>> solver;
    //Eigen::GMRES<MAT> solver;
    //Eigen::GMRES<MAT, Eigen::IdentityPreconditioner> solver;
    //Eigen::GMRES<MAT, Eigen::IncompleteLUT<double>> solver;
    //Eigen::GMRES<MAT, Eigen::DiagonalPreconditioner<double>> solver;

  public:
    DGFEMSpace1D(u_int Nx, double xl, double xr, double Nt_tol, double Nt_Ftol, double TOL);
    void BuildQuad(u_int np);
    void Projection(u_int cell, func f0, double t, bU&);
    VEC<double> Composition(const SOL&, u_int cell, double x, double t);
    void init(func f0);
    double cal_dt();
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
    EVEC NLF(const F, const SOL& sol, const SOL& soln, func,
        const double alpha, const double t, const double dt);
    void form_jacobian_rhs(const SOL& sol, const SOL& soln, const F, afunc, func,
        const double, const double, const double);
    void solve_leqn(MAT& A, const EVEC& rhs, EVEC&);
    void run(F, afunc, func, double t_end);
    int judge_positivity(const SOL&);
    void SOL2EVEC(const SOL&, EVEC&);
    void EVEC2SOL(SOL&, const EVEC&);
    VEC<double> cal_norm(const SOL&, const SOL&, int);
    VEC<double> cal_norm(const SOL& s1, int n);
    void print_solution(std::ostream&);
    MAT get_A() const;
};

#endif //DGFEMSPACE1D_H

