/*
 * File: spatial_interseptile_wmean.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 17-May-2019 14:49:50
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "spatial_interseptile_wmean.h"
#include "spatial_interseptile_wmean_emxutil.h"
#include "sort1.h"
#include "eye.h"
#include "rdivide.h"

/* Function Definitions */

/*
 * spatial_interseptile_weighted_mean: spatial filter as performed by
 *  Cartool software, instantaneous filter which removes local outliers by
 *  spatially smoothing the maps without losing its topographical
 *  characteristics: for each electrode, the values of the 6 closest
 *  neighbours are determined, plus the central electrode value itself; the 7
 *  data points are sorted; the minimal and maximal values are removed by
 *  dropping the first and last items of this list; the remaining values are
 *  then averaged, with weights proportional to the inverse distance to the
 *  central electrode. The central electrode is given a weight of 1.
 *
 *  EEG = spatial_interseptile_weighted_mean( EEG, D )
 *
 *   Inputs
 *  --------
 *  EEG: [channels x time] EEG traces
 *  D: distance matrix, e.g. result of:
 *                          >> D = squareform(pdist(XYZ));
 *
 *   Outputs
 *  ---------
 *  EEG: spatially filtered EEG
 *
 * -------------------------------------------------------------------------
 *  Cartool: https://sites.google.com/site/cartoolcommunity
 * -------------------------------------------------------------------------
 *  Renaud Marquis @ FBMlab, May 2019
 * -------------------------------------------------------------------------
 * Arguments    : emxArray_real_T *EEG
 *                emxArray_real_T *D
 * Return Type  : void
 */
void spatial_interseptile_wmean(emxArray_real_T *EEG, emxArray_real_T *D)
{
  double b_D[2];
  int i1;
  emxArray_boolean_T *r0;
  emxArray_real_T *r1;
  int ixstart;
  int trueCount;
  int ix;
  emxArray_int32_T *r2;
  int partialTrueCount;
  emxArray_real_T *Dinv;
  emxArray_int32_T *r3;
  emxArray_real_T *Neighbours;
  int c;
  emxArray_real_T *SortedDist;
  emxArray_real_T *SortedDistIdx;
  emxArray_int32_T *iidx;
  int t;
  emxArray_real_T *NewVals;
  double Cluster[7];
  double Seven[7];
  double mtmp;
  boolean_T exitg2;
  boolean_T MinSeven[7];
  boolean_T exitg1;
  int tmp_data[7];
  boolean_T MaxSeven[7];
  boolean_T b_MaxSeven;
  double M_data[7];
  int b_tmp_data[7];
  int b_trueCount;
  double Dinv_data[7];
  double x_data[7];
  double b_x_data[7];
  double y;
  int c_tmp_data[7];

  /*  you can change the number of neighbours considered (and the type of inter-xxx that will be used, e.g. interquartile, interseptile, interdecile, ...) */
  for (i1 = 0; i1 < 2; i1++) {
    b_D[i1] = D->size[i1];
  }

  emxInit_boolean_T(&r0, 2);
  emxInit_real_T(&r1, 2);
  eye(b_D, r1);
  i1 = r0->size[0] * r0->size[1];
  r0->size[0] = r1->size[0];
  r0->size[1] = r1->size[1];
  emxEnsureCapacity((emxArray__common *)r0, i1, (int)sizeof(boolean_T));
  ixstart = r1->size[0] * r1->size[1];
  for (i1 = 0; i1 < ixstart; i1++) {
    r0->data[i1] = (r1->data[i1] != 0.0);
  }

  ixstart = r0->size[0] * r0->size[1] - 1;
  trueCount = 0;
  for (ix = 0; ix <= ixstart; ix++) {
    if (r0->data[ix]) {
      trueCount++;
    }
  }

  emxInit_int32_T(&r2, 1);
  i1 = r2->size[0];
  r2->size[0] = trueCount;
  emxEnsureCapacity((emxArray__common *)r2, i1, (int)sizeof(int));
  partialTrueCount = 0;
  for (ix = 0; ix <= ixstart; ix++) {
    if (r0->data[ix]) {
      r2->data[partialTrueCount] = ix + 1;
      partialTrueCount++;
    }
  }

  ixstart = r2->size[0];
  for (i1 = 0; i1 < ixstart; i1++) {
    D->data[r2->data[i1] - 1] = rtNaN;
  }

  emxFree_int32_T(&r2);
  emxInit_real_T(&Dinv, 2);

  /*  important for both Dinv and future search of neighbours */
  rdivide(D, Dinv);

  /*  inverse distance weight matrix */
  for (i1 = 0; i1 < 2; i1++) {
    b_D[i1] = D->size[i1];
  }

  eye(b_D, r1);
  i1 = r0->size[0] * r0->size[1];
  r0->size[0] = r1->size[0];
  r0->size[1] = r1->size[1];
  emxEnsureCapacity((emxArray__common *)r0, i1, (int)sizeof(boolean_T));
  ixstart = r1->size[0] * r1->size[1];
  for (i1 = 0; i1 < ixstart; i1++) {
    r0->data[i1] = (r1->data[i1] != 0.0);
  }

  emxFree_real_T(&r1);
  ixstart = r0->size[0] * r0->size[1] - 1;
  trueCount = 0;
  for (ix = 0; ix <= ixstart; ix++) {
    if (r0->data[ix]) {
      trueCount++;
    }
  }

  emxInit_int32_T(&r3, 1);
  i1 = r3->size[0];
  r3->size[0] = trueCount;
  emxEnsureCapacity((emxArray__common *)r3, i1, (int)sizeof(int));
  partialTrueCount = 0;
  for (ix = 0; ix <= ixstart; ix++) {
    if (r0->data[ix]) {
      r3->data[partialTrueCount] = ix + 1;
      partialTrueCount++;
    }
  }

  emxFree_boolean_T(&r0);
  ixstart = r3->size[0];
  for (i1 = 0; i1 < ixstart; i1++) {
    Dinv->data[r3->data[i1] - 1] = 1.0;
  }

  emxFree_int32_T(&r3);
  emxInit_real_T(&Neighbours, 2);
  i1 = Neighbours->size[0] * Neighbours->size[1];
  Neighbours->size[0] = D->size[0];
  Neighbours->size[1] = 6;
  emxEnsureCapacity((emxArray__common *)Neighbours, i1, (int)sizeof(double));
  ixstart = D->size[0] * 6;
  for (i1 = 0; i1 < ixstart; i1++) {
    Neighbours->data[i1] = rtNaN;
  }

  c = 0;
  emxInit_real_T(&SortedDist, 2);
  emxInit_real_T(&SortedDistIdx, 2);
  emxInit_int32_T1(&iidx, 2);
  while (c <= D->size[0] - 1) {
    ixstart = D->size[1];
    i1 = SortedDist->size[0] * SortedDist->size[1];
    SortedDist->size[0] = 1;
    SortedDist->size[1] = ixstart;
    emxEnsureCapacity((emxArray__common *)SortedDist, i1, (int)sizeof(double));
    for (i1 = 0; i1 < ixstart; i1++) {
      SortedDist->data[SortedDist->size[0] * i1] = D->data[c + D->size[0] * i1];
    }

    sort(SortedDist, iidx);
    i1 = SortedDistIdx->size[0] * SortedDistIdx->size[1];
    SortedDistIdx->size[0] = 1;
    SortedDistIdx->size[1] = iidx->size[1];
    emxEnsureCapacity((emxArray__common *)SortedDistIdx, i1, (int)sizeof(double));
    ixstart = iidx->size[0] * iidx->size[1];
    for (i1 = 0; i1 < ixstart; i1++) {
      SortedDistIdx->data[i1] = iidx->data[i1];
    }

    for (i1 = 0; i1 < 6; i1++) {
      Neighbours->data[c + Neighbours->size[0] * i1] = SortedDistIdx->data[i1];
    }

    c++;
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&SortedDistIdx);
  emxFree_real_T(&SortedDist);
  i1 = EEG->size[1];
  t = 0;
  emxInit_real_T1(&NewVals, 1);
  while (t <= i1 - 1) {
    ix = NewVals->size[0];
    NewVals->size[0] = EEG->size[0];
    emxEnsureCapacity((emxArray__common *)NewVals, ix, (int)sizeof(double));
    ixstart = EEG->size[0];
    for (ix = 0; ix < ixstart; ix++) {
      NewVals->data[ix] = 0.0;
    }

    for (c = 0; c < EEG->size[0]; c++) {
      Cluster[0] = 1.0 + (double)c;
      for (ix = 0; ix < 6; ix++) {
        Cluster[ix + 1] = Neighbours->data[c + Neighbours->size[0] * ix];
      }

      for (ix = 0; ix < 7; ix++) {
        Seven[ix] = EEG->data[((int)Cluster[ix] + EEG->size[0] * t) - 1];
      }

      ixstart = 1;
      mtmp = EEG->data[((int)Cluster[0] + EEG->size[0] * t) - 1];
      if (rtIsNaN(mtmp)) {
        ix = 2;
        exitg2 = false;
        while ((!exitg2) && (ix < 8)) {
          ixstart = ix;
          if (!rtIsNaN(EEG->data[((int)Cluster[ix - 1] + EEG->size[0] * t) - 1]))
          {
            mtmp = EEG->data[((int)Cluster[ix - 1] + EEG->size[0] * t) - 1];
            exitg2 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < 7) {
        while (ixstart + 1 < 8) {
          if (EEG->data[((int)Cluster[ixstart] + EEG->size[0] * t) - 1] < mtmp)
          {
            mtmp = EEG->data[((int)Cluster[ixstart] + EEG->size[0] * t) - 1];
          }

          ixstart++;
        }
      }

      for (ix = 0; ix < 7; ix++) {
        MinSeven[ix] = (Seven[ix] == mtmp);
      }

      ixstart = 1;
      mtmp = EEG->data[((int)Cluster[0] + EEG->size[0] * t) - 1];
      if (rtIsNaN(mtmp)) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix < 8)) {
          ixstart = ix;
          if (!rtIsNaN(EEG->data[((int)Cluster[ix - 1] + EEG->size[0] * t) - 1]))
          {
            mtmp = EEG->data[((int)Cluster[ix - 1] + EEG->size[0] * t) - 1];
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < 7) {
        while (ixstart + 1 < 8) {
          if (EEG->data[((int)Cluster[ixstart] + EEG->size[0] * t) - 1] > mtmp)
          {
            mtmp = EEG->data[((int)Cluster[ixstart] + EEG->size[0] * t) - 1];
          }

          ixstart++;
        }
      }

      trueCount = 0;
      for (ix = 0; ix < 7; ix++) {
        b_MaxSeven = (Seven[ix] == mtmp);
        if (MinSeven[ix] + b_MaxSeven < 1) {
          trueCount++;
        }

        MaxSeven[ix] = b_MaxSeven;
      }

      partialTrueCount = 0;
      for (ix = 0; ix < 7; ix++) {
        if (MinSeven[ix] + MaxSeven[ix] < 1) {
          tmp_data[partialTrueCount] = ix + 1;
          partialTrueCount++;
        }
      }

      for (ix = 0; ix < trueCount; ix++) {
        M_data[ix] = EEG->data[((int)Cluster[tmp_data[ix] - 1] + EEG->size[0] *
          t) - 1];
      }

      /*  Transform a matrix into a vector using colon (:) operator */
      /*   */
      /*  (this function simply avoids using two lines of code) */
      ixstart = 0;
      for (ix = 0; ix < 7; ix++) {
        if (MinSeven[ix] + MaxSeven[ix] < 1) {
          ixstart++;
        }
      }

      partialTrueCount = 0;

      /*  Transform a matrix into a vector using colon (:) operator */
      /*   */
      /*  (this function simply avoids using two lines of code) */
      b_trueCount = 0;
      for (ix = 0; ix < 7; ix++) {
        if (MinSeven[ix] + MaxSeven[ix] < 1) {
          tmp_data[partialTrueCount] = ix + 1;
          partialTrueCount++;
        }

        if (MinSeven[ix] + MaxSeven[ix] < 1) {
          b_trueCount++;
        }
      }

      partialTrueCount = 0;
      for (ix = 0; ix < 7; ix++) {
        if (MinSeven[ix] + MaxSeven[ix] < 1) {
          b_tmp_data[partialTrueCount] = ix + 1;
          partialTrueCount++;
        }
      }

      for (ix = 0; ix < ixstart; ix++) {
        Dinv_data[ix] = Dinv->data[(int)Cluster[tmp_data[ix] - 1] - 1];
      }

      for (ix = 0; ix < trueCount; ix++) {
        x_data[ix] = M_data[ix] * Dinv_data[ix];
      }

      if (trueCount == 0) {
        mtmp = 0.0;
      } else {
        mtmp = x_data[0];
        for (ixstart = 2; ixstart <= trueCount; ixstart++) {
          mtmp += x_data[ixstart - 1];
        }
      }

      for (ix = 0; ix < b_trueCount; ix++) {
        b_x_data[ix] = Dinv->data[(int)Cluster[b_tmp_data[ix] - 1] - 1];
      }

      if (b_trueCount == 0) {
        y = 0.0;
      } else {
        for (ix = 0; ix < b_trueCount; ix++) {
          c_tmp_data[ix] = (int)Cluster[b_tmp_data[ix] - 1];
        }

        y = Dinv->data[c_tmp_data[0] - 1];
        for (ixstart = 2; ixstart <= b_trueCount; ixstart++) {
          y += b_x_data[ixstart - 1];
        }
      }

      NewVals->data[c] = mtmp / y;
    }

    ixstart = NewVals->size[0];
    for (ix = 0; ix < ixstart; ix++) {
      EEG->data[ix + EEG->size[0] * t] = NewVals->data[ix];
    }

    t++;
  }

  emxFree_real_T(&NewVals);
  emxFree_real_T(&Neighbours);
  emxFree_real_T(&Dinv);
}

/*
 * File trailer for spatial_interseptile_wmean.c
 *
 * [EOF]
 */
