/**
 * @file TemplateQuadrature.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#ifndef TEMPLATEQUADRTURE_H
#define TEMPLATEQUADRTURE_H

#include <vector>
#include <iostream>
#include <cmath>

class TemplateQuadrature {
  protected:
    u_int np;
    std::vector<double> pnt;
    std::vector<double> wei;
  public:
    std::vector<std::vector<double> > LGL(u_int np);
    void set_np(u_int np);
    void set_points(const std::vector<double>& p);
    void set_weight(const std::vector<double>& w);
    std::vector<double> points();
    std::vector<double> weight();
    void print(std::ostream&);

};

#endif

