/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: CompactSVM.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 25-Apr-2019 09:54:28
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "collision_classifier.h"
#include "CompactSVM.h"

/* Function Definitions */

/*
 * Arguments    : c_classreg_learning_coder_class *obj
 * Return Type  : void
 */
void CompactSVM_CompactSVM(c_classreg_learning_coder_class *obj)
{
  int i6;
  static const double dv15[4] = { -1.0139318885448916, 0.39512383900928133,
    49.506346749225287, 1.7558823529411685 };

  static const double dv16[4] = { 0.7097052701984341, 0.24435372788243789,
    23.859890900103604, 0.87237770699226624 };

  obj->Bias = -7.548107436699695;
  for (i6 = 0; i6 < 4; i6++) {
    obj->Mu[i6] = dv15[i6];
    obj->Sigma[i6] = dv16[i6];
  }
}

/*
 * File trailer for CompactSVM.c
 *
 * [EOF]
 */
