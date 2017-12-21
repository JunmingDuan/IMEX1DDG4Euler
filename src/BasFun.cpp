#include "BasFun.h"

std::vector<double> Poly(double x) {
  std::vector<double> poly(K);
  switch (K) {
    case 1:
      poly[0] = 1;
      return poly;
    case 2:
      poly[0] = 1;
      poly[1] = x;
      return poly;
    case 3:
      poly[0] = 1;
      poly[1] = x;
      poly[2] = (3*x*x-1)/2;
      return poly;
    case 4:
      poly[0] = 1;
      poly[1] = x;
      poly[2] = (3*x*x-1)/2;
      poly[3] = (5*x*x*x-3*x)/2;
      return poly;
    case 5:
      poly[0] = 1;
      poly[1] = x;
      poly[2] = (3*x*x-1)/2;
      poly[3] = (5*x*x*x-3*x)/2;
      poly[4] = (35*x*x*x*x-30*x*x+3)/8;
      return poly;
    default: std::cout << "Wrong basis choice!" << std::endl; abort(); break;
  }
  return poly;
}

/**
 * @brief PG Gradient of polynomial P.
 *
 * @param x
 *
 * @return A vector of length DIM.
 */
std::vector<double> PolyG(double x) {
  std::vector<double> poly(K);
  switch (K) {
    case 1:
      poly[0] = 0;
      return poly;
    case 2:
      poly[0] = 0;
      poly[1] = 1;
      return poly;
    case 3:
      poly[0] = 0;
      poly[1] = 1;
      poly[2] = 3*x;
      return poly;
    case 4:
      poly[0] = 0;
      poly[1] = 1;
      poly[2] = 3*x;
      poly[3] = (15*x*x-3)/2;
      return poly;
    case 5:
      poly[0] = 0;
      poly[1] = 1;
      poly[2] = 3*x;
      poly[3] = (15*x*x-3)/2;
      poly[4] = (140*x*x*x-60*x)/8;
      return poly;
    default: std::cout << "Wrong basis choice!" << std::endl; abort(); break;
  }
  return poly;
}

