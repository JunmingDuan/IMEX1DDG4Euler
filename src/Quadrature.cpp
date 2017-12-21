/**
 * @file Quadrature.cpp
 * @brief Quadrature points and weights on [-1,1]
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include "Quadrature.h"

void Quadrature::set_jacobi(const double j1, const double j2) {
  local_to_global_jacobian = j1;
  global_to_local_jacobian = j2;
}

double Quadrature::l2g_jacobian() {
  return local_to_global_jacobian;
}

double Quadrature::g2l_jacobian() {
  return global_to_local_jacobian;
}

void Quadrature::print(std::ostream& os) {
  os << "Number of quadrature points: " << np << std::endl;
  for(u_int i = 0; i < np; ++i) {
    os << pnt[i] << "\t";
  }
  os << "\n";
  for(u_int i = 0; i < np; ++i) {
    os << wei[i] << "\t";
  }
  os << "\n";
  os << local_to_global_jacobian << " " << global_to_local_jacobian;
  os << std::endl;
}

