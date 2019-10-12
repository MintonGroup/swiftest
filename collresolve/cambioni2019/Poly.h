/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Poly.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 25-Apr-2019 09:54:28
 */

#ifndef POLY_H
#define POLY_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "collision_classifier_types.h"

/* Function Declarations */
extern void Poly(const double svT[80], const double x[4], double kernelProduct
                 [20]);
extern void b_Poly(const double svT[184], const double x[4], double
                   kernelProduct[46]);
extern void c_Poly(const double svT[544], const double x[4], double
                   kernelProduct[136]);

#endif

/*
 * File trailer for Poly.h
 *
 * [EOF]
 */
