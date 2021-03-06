//
// File: inverseUpdatePlus.h
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 10-Aug-2017 21:27:39
//
#ifndef INVERSEUPDATEPLUS_H
#define INVERSEUPDATEPLUS_H

// Include Files
#include <cmath>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "onlineGPR_types.h"

// Function Declarations
extern void inverseUpdatePlus(const double inv_A_data[], const int inv_A_size[2],
  const double u[200], const double v[200], double inv_A_update_data[], int
  inv_A_update_size[2]);

#endif

//
// File trailer for inverseUpdatePlus.h
//
// [EOF]
//
