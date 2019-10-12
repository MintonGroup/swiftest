/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Poly.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 25-Apr-2019 09:54:28
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "collision_classifier.h"
#include "Poly.h"

/* Function Definitions */

/*
 * Arguments    : const double svT[80]
 *                const double x[4]
 *                double kernelProduct[20]
 * Return Type  : void
 */
void Poly(const double svT[80], const double x[4], double kernelProduct[20])
{
  int i0;
  double d0;
  int i1;
  for (i0 = 0; i0 < 20; i0++) {
    d0 = 0.0;
    for (i1 = 0; i1 < 4; i1++) {
      d0 += x[i1] * svT[i1 + (i0 << 2)];
    }

    kernelProduct[i0] = (d0 + 1.0) * (d0 + 1.0);
  }
}

/*
 * Arguments    : const double svT[184]
 *                const double x[4]
 *                double kernelProduct[46]
 * Return Type  : void
 */
void b_Poly(const double svT[184], const double x[4], double kernelProduct[46])
{
  int i2;
  double d1;
  int i3;
  for (i2 = 0; i2 < 46; i2++) {
    d1 = 0.0;
    for (i3 = 0; i3 < 4; i3++) {
      d1 += x[i3] * svT[i3 + (i2 << 2)];
    }

    kernelProduct[i2] = (d1 + 1.0) * (d1 + 1.0);
  }
}

/*
 * Arguments    : const double svT[544]
 *                const double x[4]
 *                double kernelProduct[136]
 * Return Type  : void
 */
void c_Poly(const double svT[544], const double x[4], double kernelProduct[136])
{
  int i4;
  double d2;
  int i5;
  for (i4 = 0; i4 < 136; i4++) {
    d2 = 0.0;
    for (i5 = 0; i5 < 4; i5++) {
      d2 += x[i5] * svT[i5 + (i4 << 2)];
    }

    kernelProduct[i4] = (d2 + 1.0) * (d2 + 1.0);
  }
}

/*
 * File trailer for Poly.c
 *
 * [EOF]
 */
