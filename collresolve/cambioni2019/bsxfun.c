/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: bsxfun.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 25-Apr-2019 09:54:28
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "collision_classifier.h"
#include "bsxfun.h"

/* Function Definitions */

/*
 * Arguments    : const double a_data[]
 *                const int a_size[2]
 *                const double b_data[]
 *                const int b_size[2]
 *                double c_data[]
 *                int c_size[2]
 * Return Type  : void
 */
void bsxfun(const double a_data[], const int a_size[2], const double b_data[],
            const int b_size[2], double c_data[], int c_size[2])
{
  int acoef;
  int bcoef;
  int sck;
  int b_b_size;
  int k;
  acoef = b_size[1];
  bcoef = a_size[1];
  if (acoef < bcoef) {
    bcoef = acoef;
  }

  if (b_size[1] == 1) {
    sck = a_size[1];
  } else if (a_size[1] == 1) {
    sck = b_size[1];
  } else if (a_size[1] == b_size[1]) {
    sck = a_size[1];
  } else {
    sck = bcoef;
  }

  c_size[0] = 1;
  acoef = b_size[1];
  bcoef = a_size[1];
  if (acoef < bcoef) {
    bcoef = acoef;
  }

  if (b_size[1] == 1) {
    c_size[1] = (signed char)a_size[1];
  } else if (a_size[1] == 1) {
    c_size[1] = (signed char)b_size[1];
  } else if (a_size[1] == b_size[1]) {
    c_size[1] = (signed char)a_size[1];
  } else {
    c_size[1] = (signed char)bcoef;
  }

  acoef = b_size[1];
  bcoef = a_size[1];
  if (acoef < bcoef) {
    bcoef = acoef;
  }

  if (b_size[1] == 1) {
    b_b_size = a_size[1];
  } else if (a_size[1] == 1) {
    b_b_size = b_size[1];
  } else if (a_size[1] == b_size[1]) {
    b_b_size = a_size[1];
  } else {
    b_b_size = bcoef;
  }

  if ((signed char)b_b_size != 0) {
    acoef = (a_size[1] != 1);
    bcoef = (b_size[1] != 1);
    for (k = 0; k < (signed char)sck; k++) {
      c_data[k] = a_data[a_size[0] * (acoef * k)] / b_data[b_size[0] * (bcoef *
        k)];
    }
  }
}

/*
 * File trailer for bsxfun.c
 *
 * [EOF]
 */
