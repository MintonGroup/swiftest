/**
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 25-Apr-2019 09:54:28
 *
 * @file
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "collision_classifier.h"
#include "minOrMax.h"
#include "CompactClassificationECOC.h"

/* Function Definitions */

/**
 * Classifier of collision outcome for a giant impact.
 *
 * @param X predictors. Array of 4 values:
 *          - log10( mass of the target / Earth mass )
 *          - ratio mass projectile / mass target
 *          - impact angle [degree]
 *          - ration impact velocity / escape velocity
 * @param label Collision label:
 *          - 1 => Hit-and-run
 *          - 0 => Accretion
 *          - -1 => Disruption
 * @param score SVM scores for each class
 */
void collision_classifier(const double X[4], int* label, double score[3]) {
	int idx;
	signed char CompactMdl_CodingMatrix[9];
	static const signed char iv0[9] = { 1, -1, 0, 1, 0, -1, 0, 1, -1 };

	double CompactMdl_Prior[3];
	static const double dv0[3] = { 0.091420534458509145, 0.37130801687763715, 0.53727144866385368 };

	double pscore[3];
	double classnum;
	boolean_T notNaN;
	double y;
	double M[9];
	int c;
	int k;

	for (idx = 0; idx < 9; idx++) {
		CompactMdl_CodingMatrix[idx] = iv0[idx];
	}

	for (idx = 0; idx < 3; idx++) {
		CompactMdl_Prior[idx] = dv0[idx];
	}

	for (idx = 0; idx < 9; idx++) {
		classnum = CompactMdl_CodingMatrix[idx];
		if (CompactMdl_CodingMatrix[idx] == 0) {
			classnum = rtNaN;
		}

		M[idx] = classnum;
	}

	notNaN = false;
	localScore(X, pscore);
	for (idx = 0; idx < 3; idx++) {
		y = 0.0;
		c = 0;
		for (k = 0; k < 3; k++) {
			classnum = 1.0 - M[idx + 3 * k] * pscore[k];
			if ((0.0 > classnum) || rtIsNaN(classnum)) {
				classnum = 0.0;
			}

			if (!rtIsNaN(classnum)) {
				y += classnum;
				c++;
			}
		}

		if (c == 0) {
			y = rtNaN;
		} else {
			y /= (double)c;
		}

		classnum = -(y / 2.0);
		score[idx] = classnum;

		if (!rtIsNaN(classnum)) {
			notNaN = true;
		}
	}

	eml_extremum(CompactMdl_Prior, &classnum, &idx);
	classnum = rtNaN;
	if (notNaN) {
		eml_extremum(score, &classnum, &c);
		classnum = c;

		if (classnum < 4.294967296E+9) {
			if (classnum >= 0.0) {
				idx = (int)classnum;
			} else {
				idx = 0;
			}
		} else {
			idx = 0;
		}
	}

	*label = idx - 2;
}

/*
 * File trailer for collision_classifier.c
 *
 * [EOF]
 */
