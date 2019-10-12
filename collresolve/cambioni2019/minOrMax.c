/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: minOrMax.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 25-Apr-2019 09:54:28
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "collision_classifier.h"
#include "minOrMax.h"

/* Function Definitions */

/*
 * Arguments    : const double x[3]
 *                double *extremum
 *                int *indx
 * Return Type  : void
 */
void eml_extremum(const double x[3], double *extremum, int *indx)
{
  int ixstart;
  int ix;
  boolean_T exitg1;
  ixstart = 1;
  *extremum = x[0];
  *indx = 1;
  if (rtIsNaN(x[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 4)) {
      ixstart = ix;
      if (!rtIsNaN(x[ix - 1])) {
        *extremum = x[ix - 1];
        *indx = ix;
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 3) {
    while (ixstart + 1 < 4) {
      if (x[ixstart] > *extremum) {
        *extremum = x[ixstart];
        *indx = ixstart + 1;
      }

      ixstart++;
    }
  }
}

/*
 * File trailer for minOrMax.c
 *
 * [EOF]
 */
