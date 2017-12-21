#ifndef INTERVAL_CRD_TRS_H
#define INTERVAL_CRD_TRS_H
/**
 * @file interval_crd_trs.h
 * @brief 1D coordinate transformation
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-04
 */

void local_to_global(const double lp, const std::vector<double>& lv, const std::vector<double>& gv, double* gp); //local points -> global points.

void global_to_local(const double gp, const std::vector<double>& lv, const std::vector<double>& gv, double* lp); //global points -> local points.

double local_to_global_jacobian(const std::vector<double>& lv, const std::vector<double>& gv);

double global_to_local_jacobian(const std::vector<double>& lv, const std::vector<double>& gv);

#endif //INTERVAL_CRD_TRS_H

