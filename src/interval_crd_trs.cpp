#include <vector>

void local_to_global(const double lp, const std::vector<double>& lv, const std::vector<double>& gv, double* gp) //local points -> global points.
{
  double lambda[2];
  lambda[0] = (lv[1] - lp)/(lv[1] - lv[0]);
  lambda[1] = (lp - lv[0])/(lv[1] - lv[0]);
  *gp = lambda[0]*gv[0] + lambda[1]*gv[1];
}

void global_to_local(const double gp, const std::vector<double>& lv, const std::vector<double>& gv, double* lp) //global points -> local points.
{
	double lambda[2];
	lambda[0] = (gv[1] - gp)/(gv[1] - gv[0]);
	lambda[1] = (gp - gv[0])/(gv[1] - gv[0]);
	*lp = lambda[0]*lv[0] + lambda[1]*lv[1];
}

double local_to_global_jacobian(const std::vector<double>& lv, const std::vector<double>& gv) //local -> global jacobian.
{
	return (gv[1] - gv[0])/(lv[1] - lv[0]);
}

double global_to_local_jacobian(const std::vector<double>& lv, const std::vector<double>& gv) //global -> local jacobian.
{
	return (lv[1] - lv[0])/(gv[1] - gv[0]);
}

