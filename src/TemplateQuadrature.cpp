/**
 * @file Quadrature.cpp
 * @brief Quadrature points and weights on [-1,1]
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include "TemplateQuadrature.h"

std::vector<std::vector<double> > TemplateQuadrature::LGL(u_int np) {
  std::vector<std::vector<double> > pw;
  pw.resize(2);
  pw[0].resize(np);
  pw[1].resize(np);
  switch (np) { //from wiki
    case 2: {
              pw[0][0] = -1, pw[0][1] = 1;
              pw[1][0] = 1, pw[1][1] = 1;
              break;
            }
    case 3: {
              pw[0][0] = -1, pw[0][1] = 0, pw[0][2] = 1;
              pw[1][0] = 1./3, pw[1][1] = 4./3, pw[1][2] = 1./3;
              break;
            }
    case 4: {
              pw[0][0] = -1, pw[0][1] = -1/sqrt(5), pw[0][2] = 1/sqrt(5), pw[0][3] = 1;
              pw[1][0] = 1./6,  pw[1][1] = 5./6, pw[1][2] = 5./6, pw[1][3] = 1./6;
              break;
            }
    case 5: {
              pw[0][0] = -1, pw[0][1] = -sqrt(3./7), pw[0][2] = 0, pw[0][3] = sqrt(3./7), pw[0][4] = 1;
              pw[1][0] = 1./10, pw[1][1] = 49./90, pw[1][2] = 32./45, pw[1][3] = 49./90, pw[1][4] = 1./10;
              break;
            }
    case 6: {
              pw[0][0] = -1, pw[0][1] = -sqrt(1./3+2*sqrt(7)/21), pw[0][2] = -sqrt(1./3-2*sqrt(7)/21),
              pw[0][3] = sqrt(1./3-2*sqrt(7)/21), pw[0][4] = sqrt(1./3+2*sqrt(7)/21), pw[0][5] = 1;
              pw[1][0] = 1./15, pw[1][1] = (14-sqrt(7))/30, pw[1][2] = (14+sqrt(7))/30,
              pw[1][3] = (14+sqrt(7))/30, pw[1][4] = (14-sqrt(7))/30, pw[1][5] = 1./15;
              break;
            }
    case 7: {
              pw[0][0] = -1, pw[0][1] = -sqrt(5./11+2./11*sqrt(5./3)), pw[0][2] = -sqrt(5./11-2./11*sqrt(5./3)),
              pw[0][3] = 0,
              pw[0][4] = sqrt(5./11-2./11*sqrt(5./3)), pw[0][5] = sqrt(5./11+2./11*sqrt(5./3)), pw[0][6] = 1;
              pw[1][0] = 1./21, pw[1][1] = (124-7*sqrt(15))/350, pw[1][2] = (124+7*sqrt(15))/350,
              pw[1][3] = 256./525,
              pw[1][4] = (124+7*sqrt(15))/350, pw[1][5] = (124-7*sqrt(15))/350, pw[1][6] = 1./21;
              break;
            }
    default:{
               std::cout << "Wrong choice!" << std::endl;
               break;
            }
  }
  return pw;
}

void TemplateQuadrature::set_np(u_int n) {
  np = n;
}

void TemplateQuadrature::set_weight(const std::vector<double>& w) {
  wei = w;
}

void TemplateQuadrature::set_points(const std::vector<double>& p) {
  pnt = p;
}

std::vector<double> TemplateQuadrature::points() {
  return pnt;
}

std::vector<double> TemplateQuadrature::weight() {
  return wei;
}

void TemplateQuadrature::print(std::ostream& os) {
  os << "Number of quadrature points: " << np << std::endl;
  for(u_int i = 0; i < np; ++i) {
    os << pnt[i] << "\t";
  }
  os << "\n";
  for(u_int i = 0; i < np; ++i) {
    os << wei[i] << "\t";
  }
  os << std::endl;
}

