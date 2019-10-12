/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: collision_classifier_types.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 25-Apr-2019 09:54:28
 */

#ifndef COLLISION_CLASSIFIER_TYPES_H
#define COLLISION_CLASSIFIER_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_c_classreg_learning_coder_class
#define typedef_c_classreg_learning_coder_class

typedef struct {
  double Alpha[136];
  double Bias;
  double Mu[4];
  double Sigma[4];
  double ClassNames[2];
  int ClassNamesLength[2];
  double Prior[2];
  double NonzeroProbClasses[2];
  double Cost[4];
} c_classreg_learning_coder_class;

#endif                                 /*typedef_c_classreg_learning_coder_class*/
#endif

/*
 * File trailer for collision_classifier_types.h
 *
 * [EOF]
 */
