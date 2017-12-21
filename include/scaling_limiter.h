#ifndef SCALING_LIMITER_H
#define SCALING_LIMITER_H

#include <algorithm>
#include "vvector.h"
#include "para.h"

#define VECD VEC<double>

class scaling_limiter {
  private:
    //VEC<VECD> cell_average;
    //VEC<VEC<VECD>> val_int_pnt;
  public:
    //scaling_limiter();
    //scaling_limiter(const VEC<VECD>& av, VEC<VEC<VECD>>& val);
    //void set(const VEC<VECD>& av, VEC<VEC<VECD>>& val);
    static void run(const VEC<VECD>& av, const VEC<VEC<VECD>>& val, SOL& sol);
};

#endif //SCALING_LIMITER_H

