/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Reconstruct an array along a given direction.

  General purpose reconstruction routine using limited linear
  interpolation.
  Given a 3D array q, it reconstruct along the "dir" direction.

  \authors A. Mignone (mignone@to.infn.it)\n

  \date    Jul 19, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void ArrayReconstruct(double ***Q, char *flag, int i, int j, int k, int dir,
                      double *qL, double *qR, int rec, Grid *grid)
/*!
 *********************************************************************** */
{
  double dqp, dqm, dq;
  PLM_Coeffs plm_coeffs;
  int rec0 = (rec/10)*10;
  int ntot;
  static double *qfwd, *qbck;

  if (qfwd == NULL){
    qfwd = ARRAY_1D(NMAX_POINT, double);
    qbck = ARRAY_1D(NMAX_POINT, double); 
  }

  if (dir == IDIR) {
    ntot = NX1_TOT;
    for (i = 0; i < ntot; i++) qfwd[i] = Q[k][j][i];
  } 
 
  if (dir == JDIR) {
    ntot = NX2_TOT;
    for (j = 0; j < ntot; j++) qfwd[j] = Q[k][j][i];
  }

  if (dir == KDIR) {
    ntot = NX3_TOT;
    for (k = 0; k < ntot; k++) qfwd[k] = Q[k][j][i];
  }

/* ----------------------------------------------
   Flat reconstruction
   ---------------------------------------------- */
  
  if (rec0 == FLAT) {
    int n;
    dqm = qfwd[1] - qfwd[0];
    for (n = 1; n < ntot-1; n++){      
      qL[n]   = qfwd[n];
      qR[n-1] = qfwd[n];
    }
    return;
  }

/* ----------------------------------------------
   Linear reconstruction
   ---------------------------------------------- */
  
  if (rec0 == LINEAR) {
    int n;
    dqm = qfwd[1] - qfwd[0];
    for (n = 1; n < ntot-1; n++){
      
      dqp = qfwd[n+1] - qfwd[n];
      dqm = qfwd[n] - qfwd[n-1];

      if      (rec == FLAT_LIM)    dq = 0.0;
      else if (rec == MINMOD_LIM)  dq = MINMOD_LIMITER(dqp, dqm);
      else if (rec == VANLEER_LIM) dq = VANLEER_LIMITER(dqp, dqm);
      else if (rec == MC_LIM)      dq = MC_LIMITER(dqp, dqm);
      else{
        SET_LIMITER (dq, dqp, dqm, 2.0, 2.0);
      }
      qL[n]   = qfwd[n] + 0.5*dq;
      qR[n-1] = qfwd[n] - 0.5*dq;
      
    }
    return;
  }

/* ----------------------------------------------
   PARABOLIC reconstruction
   ---------------------------------------------- */

  if (rec == PARABOLIC) {
    for (i = 0; i < ntot; i++) qbck[i] = qfwd[ntot-i-1];

    int beg = 1, end = ntot-2;
    for (i = beg; i <= end; i++){
 
      double *v = qfwd;
      double dv, dvp, dvm;
      double dvc, Sm1, Sp1, Sp2, SM;
      double ap, am, vp, vm;

      #if SHOCK_FLATTENING == MULTID
      if (flag[i] & FLAG_MINMOD) {
        dqp = v[i+1] - v[i];
        dqm = v[i] - v[i-1];
        dq = MINMOD_LIMITER(dqp, dqm);
        qL[i]   = v[i] + 0.5*dq;
        qR[i-1] = v[i] - 0.5*dq;
        continue;
      }
      #endif
  
      vp = (    -v[i-1] + 5.0*v[i] + 2.0*v[i+1])/6.0;
      vm = ( 2.0*v[i-1] + 5.0*v[i] -     v[i+1])/6.0;
      
      dvp = vp - v[i];
      dvm = vm - v[i];
  
      dv  = v[i+1] - v[i];
      vp  = v[i] + MINMOD_LIMITER(dvp, dv);
      
      dv  = v[i] - v[i-1];
      vm  = v[i] + MINMOD_LIMITER(dvm, -dv);
       
      dvp = vp - v[i];
      dvm = vm - v[i];
  
      if (dvp*dvm >= 0.0) dvp = dvm = 0.0;
      else{
        if      (fabs(dvp) >= 2.0*fabs(dvm)) dvp = -2.0*dvm;
        else if (fabs(dvm) >= 2.0*fabs(dvp)) dvm = -2.0*dvp;
      }
      qL[i]   = vp = v[i] + dvp; 
      qR[i-1] = vm = v[i] + dvm;
  
    }
    return;
  }

/* ----------------------------------------------
   WENO3 reconstruction
   ---------------------------------------------- */

  if (rec == WENO3) {

    for (i = 0; i < ntot; i++) qbck[i] = qfwd[ntot-i-1];
  
    int beg = 1, end = ntot-2;
    for (i = beg; i <= end; i++){
      int ip = i;
      int im = end-(i-beg);
      #if SHOCK_FLATTENING == MULTID
      if (flag[i] & FLAG_MINMOD) {
        dqp = qfwd[i+1] - qfwd[i];
        dqm = qfwd[i] - qfwd[i-1];
        dq = MINMOD_LIMITER(dqp, dqm);
        qL[i]   = qfwd[i] + 0.5*dq;
        qR[i-1] = qfwd[i] - 0.5*dq;
        continue;
      }
      #endif
  
      qL[i]   = WENO3_Reconstruct (qfwd, grid->dx[dir][i], ip);
      qR[i-1] = WENO3_Reconstruct (qbck, grid->dx[dir][i], im);
    }
    return;
  }

/* ----------------------------------------------
   WENOZ reconstruction
   ---------------------------------------------- */

#if RECONSTRUCTION == WENOZ
  if (rec == WENOZ) {

    int beg = 2, end = ntot-3;
    for (i = beg; i <= end; i++){
      #if SHOCK_FLATTENING == MULTID
      if (flag[i] & FLAG_MINMOD) {
        dqp = qfwd[i+1] - qfwd[i];
        dqm = qfwd[i] - qfwd[i-1];
        dq = MINMOD_LIMITER(dqp, dqm);
        qL[i]   = qfwd[i] + 0.5*dq;
        qR[i-1] = qfwd[i] - 0.5*dq;
        continue;
      }
      #endif

      qL[i]   = WENOZ_States (qfwd, i, +1);
      qR[i-1] = WENOZ_States (qfwd, i, -1);
    }
    return;
  }
#endif
  
/* ----------------------------------------------
   MP5 reconstruction
   ---------------------------------------------- */

#if RECONSTRUCTION == MP5
  if (rec == MP5) {

    int beg = 2, end = ntot-3;
    for (i = beg; i <= end; i++){
      #if SHOCK_FLATTENING == MULTID
      if (flag[i] & FLAG_MINMOD) {
        dqp = qfwd[i+1] - qfwd[i];
        dqm = qfwd[i] - qfwd[i-1];
        dq = MINMOD_LIMITER(dqp, dqm);
        qL[i]   = qfwd[i] + 0.5*dq;
        qR[i-1] = qfwd[i] - 0.5*dq;
        continue;
      }
      #endif

      qL[i]   = MP5_States (qfwd, i, +1);
      qR[i-1] = MP5_States (qfwd, i, -1);
    }
    return;
  }
#endif
  
  printLog ("! ArrayReconstruct(): invalid rconstruction type (rec = %d)\n",
               rec);
  QUIT_PLUTO(1);

}
