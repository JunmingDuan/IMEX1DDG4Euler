/**
 * @file Quadrature.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#ifndef QUADRTURE_H
#define QUADRTURE_H

#include "TemplateQuadrature.h"

class Quadrature : public TemplateQuadrature {
  private:
    double local_to_global_jacobian;
    double global_to_local_jacobian;
  public:
    void set_jacobi(const double j1, const double j2);
    double l2g_jacobian();
    double g2l_jacobian();
    void print(std::ostream&);
};

#endif //QUADRTURE_H

