/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Collects different EMF averaging schemes.

  \author A. Mignone
  \date   March 20, 2020
*/  
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CT_EMF_ArithmeticAverage (const EMF *Z1, const double w)
/*!
 *  Compute arithmetic average of EMF at cell edges by combining
 *  the four electric field values computed at zone faces 
 *  as upwind Godunov fluxes into an edge-centered value.\n
 *  The face-centered EMF should have been stored by previous calls
 *  to CT_StoreEMF() during the one-dimensional sweeps.\n
 *  This function employs a simple arithmetic averaging of the 
 *  face-centered electric field.
 *
 * \b References:
 *    - "A Staggered Mesh Algorithm Using High Order Godunov Fluxes to 
 *       Ensure Solenoidal Magnetic Fields in Magnetohydrodynamic 
 *       Simulations"\n
 *       Balsara \& Spicer, JCP (1999) 149, 270.
 *
 * \param [in]      Z1    pointer to EMF structure
 * \param [in]      w     weighting factor
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int i, j, k;

  for (k = Z1->kbeg; k <= Z1->kend; k++){
  for (j = Z1->jbeg; j <= Z1->jend; j++){
  for (i = Z1->ibeg; i <= Z1->iend; i++){      
    #if INCLUDE_JDIR && INCLUDE_KDIR
    Z1->Ex1e[k][j][i] = w*(  Z1->exk[k][j][i] + Z1->exk[k][j+1][i] 
                           + Z1->exj[k][j][i] + Z1->exj[k+1][j][i]);
    #endif

    #if INCLUDE_IDIR && INCLUDE_KDIR
    Z1->Ex2e[k][j][i] = w*(  Z1->eyi[k][j][i] + Z1->eyi[k+1][j][i] 
                           + Z1->eyk[k][j][i] + Z1->eyk[k][j][i+1]);
    #endif

    #if INCLUDE_IDIR && INCLUDE_JDIR
    Z1->Ex3e[k][j][i] = w*(  Z1->ezi[k][j][i] + Z1->ezi[k][j+1][i] 
                           + Z1->ezj[k][j][i] + Z1->ezj[k][j][i+1]);
    #endif

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    #if DIMENSIONS == 3
    Z1->Bx1e[k][j][i] = w*(  Z1->Bxk[k][j][i] + Z1->Bxk[k][j+1][i] 
                           + Z1->Bxj[k][j][i] + Z1->Bxj[k+1][j][i]);
    Z1->Bx2e[k][j][i] = w*(  Z1->Byi[k][j][i] + Z1->Byi[k+1][j][i] 
                           + Z1->Byk[k][j][i] + Z1->Byk[k][j][i+1]);
    #endif 
    Z1->Bx3e[k][j][i] = w*(  Z1->Bzi[k][j][i] + Z1->Bzi[k][j+1][i] 
                           + Z1->Bzj[k][j][i] + Z1->Bzj[k][j][i+1]);
    #endif
    
  }}}
}

/* ********************************************************************* */
void CT_EMF_IntegrateToCorner (const Data *d, const EMF *emf, Grid *grid)
/*!
 *  Add derivatives to the 4-point arithmetic average of magnetic fields.
 *  Obtain the electric field at corners. 
 *
 * \b References: 
 *  - "An unsplit Godunov method for ideal MHD via constrained transport"\n
 *    Gardiner & Stone, JCP (2005) 205, 509.
 *    See Eq. (41), (45) and (50).
 *
 * \param [in]      d     pointer to PLUTO Data structure
 * \param [in]      emf   pointer to EMF structure
 * \param [in]      grid  pointer to Grid structure
 *
 * \return  This function has no return value.
 *********************************************************************** */
#define DEX_DYP(k,j,i) (emf->exj[k][j][i] - d->Ex1[k][j][i])
#define DEX_DZP(k,j,i) (emf->exk[k][j][i] - d->Ex1[k][j][i])

#define DEY_DXP(k,j,i) (emf->eyi[k][j][i] - d->Ex2[k][j][i])
#define DEY_DZP(k,j,i) (emf->eyk[k][j][i] - d->Ex2[k][j][i])

#define DEZ_DXP(k,j,i) (emf->ezi[k][j][i] - d->Ex3[k][j][i])
#define DEZ_DYP(k,j,i) (emf->ezj[k][j][i] - d->Ex3[k][j][i])

#define DEX_DYM(k,j,i) (d->Ex1[k][j][i] - emf->exj[k][j-1][i])
#define DEX_DZM(k,j,i) (d->Ex1[k][j][i] - emf->exk[k-1][j][i])

#define DEY_DXM(k,j,i) (d->Ex2[k][j][i] - emf->eyi[k][j][i-1])
#define DEY_DZM(k,j,i) (d->Ex2[k][j][i] - emf->eyk[k-1][j][i])

#define DEZ_DXM(k,j,i) (d->Ex3[k][j][i] - emf->ezi[k][j][i-1])
#define DEZ_DYM(k,j,i) (d->Ex3[k][j][i] - emf->ezj[k][j-1][i])
{
  int    i, j, k;
  int    iu, ju, ku;
  signed char  sx, sy, sz;

#ifdef CTU
  if (g_intStage == 1) return; /* -- not needed in predictor step of CTU -- */
#endif

  for (k = emf->kbeg; k <= emf->kend + INCLUDE_KDIR; k++){
  for (j = emf->jbeg; j <= emf->jend + INCLUDE_JDIR; j++){
  for (i = emf->ibeg; i <= emf->iend + INCLUDE_IDIR; i++){      

    DIM_EXPAND(sx = emf->svx[k][j][i];  ,
             sy = emf->svy[k][j][i];  ,
             sz = emf->svz[k][j][i];)

    DIM_EXPAND(iu = sx > 0 ? i:i+1;  ,  /* -- upwind index -- */
             ju = sy > 0 ? j:j+1;  ,
             ku = sz > 0 ? k:k+1;)

 /* ----------------------------------------
      Span X - Faces:    dEz/dy, dEy/dz
    ---------------------------------------- */

    if (sx == 0) {
      #if INCLUDE_IDIR && INCLUDE_JDIR
      emf->Ex3e[k][j][i]   += 0.5*(DEZ_DYP(k,j,i) + DEZ_DYP(k,j,i+1));
      emf->Ex3e[k][j-1][i] -= 0.5*(DEZ_DYM(k,j,i) + DEZ_DYM(k,j,i+1));
      #endif
      #if INCLUDE_IDIR && INCLUDE_KDIR
      emf->Ex2e[k][j][i]   += 0.5*(DEY_DZP(k,j,i) + DEY_DZP(k,j,i+1));
      emf->Ex2e[k-1][j][i] -= 0.5*(DEY_DZM(k,j,i) + DEY_DZM(k,j,i+1));
      #endif
    }else{
      #if INCLUDE_IDIR && INCLUDE_JDIR
      emf->Ex3e[k][j][i]   += DEZ_DYP(k,j,iu);
      emf->Ex3e[k][j-1][i] -= DEZ_DYM(k,j,iu);
      #endif
      #if INCLUDE_IDIR && INCLUDE_KDIR
      emf->Ex2e[k][j][i]   += DEY_DZP(k,j,iu);
      emf->Ex2e[k-1][j][i] -= DEY_DZM(k,j,iu);
      #endif
    }

 /* ----------------------------------------
      Span Y - Faces:    dEz/dx, dEx/dz
    ---------------------------------------- */

    if (sy == 0) {
      #if INCLUDE_IDIR && INCLUDE_JDIR
      emf->Ex3e[k][j][i]   += 0.5*(DEZ_DXP(k,j,i) + DEZ_DXP(k,j+1,i));
      emf->Ex3e[k][j][i-1] -= 0.5*(DEZ_DXM(k,j,i) + DEZ_DXM(k,j+1,i));
      #endif
      #if INCLUDE_JDIR && INCLUDE_KDIR
      emf->Ex1e[k][j][i]   += 0.5*(DEX_DZP(k,j,i) + DEX_DZP(k,j+1,i));
      emf->Ex1e[k-1][j][i] -= 0.5*(DEX_DZM(k,j,i) + DEX_DZM(k,j+1,i));
      #endif
    }else{
      #if INCLUDE_IDIR && INCLUDE_JDIR
      emf->Ex3e[k][j][i]   += DEZ_DXP(k,ju,i);
      emf->Ex3e[k][j][i-1] -= DEZ_DXM(k,ju,i);
      #endif
      #if INCLUDE_JDIR && INCLUDE_KDIR
      emf->Ex1e[k][j][i]   += DEX_DZP(k,ju,i);
      emf->Ex1e[k-1][j][i] -= DEX_DZM(k,ju,i);
      #endif
    }

 /* ----------------------------------------
      Span Z - Faces:    dEx/dy, dEy/dx
    ---------------------------------------- */

#if DIMENSIONS == 3
    if (sz == 0) {
      #if INCLUDE_JDIR && INCLUDE_KDIR
      emf->Ex1e[k][j][i]   += 0.5*(DEX_DYP(k,j,i) + DEX_DYP(k+1,j,i));
      emf->Ex1e[k][j-1][i] -= 0.5*(DEX_DYM(k,j,i) + DEX_DYM(k+1,j,i));
      #endif
      #if INCLUDE_IDIR && INCLUDE_KDIR
      emf->Ex2e[k][j][i]   += 0.5*(DEY_DXP(k,j,i) + DEY_DXP(k+1,j,i));
      emf->Ex2e[k][j][i-1] -= 0.5*(DEY_DXM(k,j,i) + DEY_DXM(k+1,j,i));
      #endif
    }else{
      #if INCLUDE_JDIR && INCLUDE_KDIR
      emf->Ex1e[k][j][i]   += DEX_DYP(ku,j,i);
      emf->Ex1e[k][j-1][i] -= DEX_DYM(ku,j,i);
      #endif
      #if INCLUDE_IDIR && INCLUDE_KDIR
      emf->Ex2e[k][j][i]   += DEY_DXP(ku,j,i);
      emf->Ex2e[k][j][i-1] -= DEY_DXM(ku,j,i);
      #endif
    }
#endif
  }}}
}
#undef DEX_DYP
#undef DEX_DZP
#undef DEY_DXP
#undef DEY_DZP
#undef DEZ_DXP
#undef DEZ_DYP
#undef DEX_DYM
#undef DEX_DZM
#undef DEY_DXM
#undef DEY_DZM
#undef DEZ_DXM
#undef DEZ_DYM

/* ********************************************************************* */
void CT_EMF_Riemann2D(const Data *d, const EMF *emf, Grid *grid)
/*!
 * Reconstruct the electric field to zone edge using the general
 * UCT formalism,
 *
 * E = (aL*vL[i]*bL[i] + aR*vR[i]*bR[i]) +- (dR*bR - dL*bL)
 *
 * where the a and d coefficients are retrieved from the 1D Riemann solver
 * at zone edges.
 * This function is used by UCT_HLL, UCT_HLLD and UCT_GFORCE
 *
 *********************************************************************** */
#define UPWIND_AVERAGE  NO
{
  int i, j, k;
#if RECONSTRUCTION == LINEAR
  #if LIMITER == DEFAULT
  int recV = VANLEER_LIM;
  #else
  int recV = LIMITER;
  #endif
  int recB = MC_LIM;
#else
  int recV = RECONSTRUCTION;
  int recB = RECONSTRUCTION;
#endif
  DIM_EXPAND(double ***Bx1s = d->Vs[BX1s];  ,
             double ***Bx2s = d->Vs[BX2s];  ,
             double ***Bx3s = d->Vs[BX3s];)
#if BACKGROUND_FIELD == YES
  double B0[3];
  double *x  = grid->x[IDIR],  *y  = grid->x[JDIR],  *z  = grid->x[KDIR];
  double *xr = grid->xr[IDIR], *yr = grid->xr[JDIR], *zr = grid->xr[KDIR];
#endif

  static char   *flag;
  static double *vL, *vR, *bL, *bR;
  double SL, SR, aL, aR, dL, dR, phi;
#if UPWIND_AVERAGE == YES
  double alphaR, alphaL, wN, wS, scrh;
#endif

  if (bL == NULL) {
    vL = ARRAY_1D(NMAX_POINT, double);
    vR = ARRAY_1D(NMAX_POINT, double);
    bL = ARRAY_1D(NMAX_POINT, double);
    bR = ARRAY_1D(NMAX_POINT, double);
    flag = ARRAY_1D(NMAX_POINT, char);
  }

/* --------------------------------------------------------
   X1. Reconstruct along the x1-direction:

     o -vx(j+1/2) and By(j+1/2) --> (i+1/2, j+1/2, k)
     o +vx(k+1/2) and Bz(k+1/2) --> (i+1/2, j, k+1/2)
   -------------------------------------------------------- */
   
  #if INCLUDE_IDIR
  for (k = emf->kbeg; k <= emf->kend; k++){ 
  for (j = emf->jbeg; j <= emf->jend; j++){

    #if INCLUDE_JDIR
    #if SHOCK_FLATTENING == MULTID
    for (i = emf->ibeg; i <= emf->iend+1; i++){
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k][j+1][i] & FLAG_MINMOD) ) flag[i] |= FLAG_MINMOD;
      else flag[i] = 0;
    }
    #endif
    ArrayReconstruct (emf->ezj, flag, i, j, k, IDIR, vL, vR, recV, grid); /* -vx */
    ArrayReconstruct (Bx2s,     flag, i, j, k, IDIR, bL, bR, recB, grid); 
    for (i = emf->ibeg; i <= emf->iend; i++){
      #if BACKGROUND_FIELD == YES
      BackgroundField (xr[i], yr[j], z[k], B0);
      bR[i] += B0[JDIR];
      bL[i] += B0[JDIR];
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      aL = aR = 0.5;
      dL = 0.5*(emf->dxL[k][j][i] + emf->dxL[k][j+1][i]);
      dR = 0.5*(emf->dxR[k][j][i] + emf->dxR[k][j+1][i]);
      #else
      #if UPWIND_AVERAGE == YES
      alphaR = MAX(emf->SyR[k][j][i], emf->SyR[k][j][i+1]);
      alphaR = MAX(0,alphaR);
      alphaL = MIN(emf->SyL[k][j][i], emf->SyL[k][j][i+1]);
      alphaL = MIN(0,alphaL);
      
      scrh = 1.0/(alphaR - alphaL);
      wS =  alphaR*scrh;
      wN = -alphaL*scrh;
      aL = wS*emf->axL[k][j][i] + wN*emf->axL[k][j+1][i];
      aR = wS*emf->axR[k][j][i] + wN*emf->axR[k][j+1][i];
      dL = wS*emf->dxL[k][j][i] + wN*emf->dxL[k][j+1][i];
      dR = wS*emf->dxR[k][j][i] + wN*emf->dxR[k][j+1][i];
      #else
      aL = 0.5*(emf->axL[k][j][i] + emf->axL[k][j+1][i]);
      aR = 0.5*(emf->axR[k][j][i] + emf->axR[k][j+1][i]);
      dL = 0.5*(emf->dxL[k][j][i] + emf->dxL[k][j+1][i]);
      dR = 0.5*(emf->dxR[k][j][i] + emf->dxR[k][j+1][i]);
      #endif
      #endif

      phi = dR*bR[i] - dL*bL[i];
      emf->Ex3e[k][j][i] = (aL*vL[i]*bL[i] + aR*vR[i]*bR[i]) + phi;

    }
    #endif  /* INCLUDE_JDIR */

    #if INCLUDE_KDIR
    #if SHOCK_FLATTENING == MULTID
    for (i = emf->ibeg; i <= emf->iend+1; i++){
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k+1][j][i] & FLAG_MINMOD) ) flag[i] |= FLAG_MINMOD;
      else flag[i] = 0;
    }
    #endif
    ArrayReconstruct (emf->eyk, flag, i, j, k, IDIR, vL, vR, recV, grid); /* +vx */
    ArrayReconstruct (Bx3s,     flag, i, j, k, IDIR, bL, bR, recB, grid); 
    for (i = emf->ibeg; i <= emf->iend; i++){
      #if BACKGROUND_FIELD == YES
      BackgroundField (xr[i], yr[j], z[k], B0);
      bR[i] += B0[KDIR];
      bL[i] += B0[KDIR];
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      aL  = aR = 0.5;
      dL = 0.5*(emf->dxL[k][j][i] + emf->dxL[k+1][j][i]);
      dR = 0.5*(emf->dxR[k][j][i] + emf->dxR[k+1][j][i]);
      #else
      #if UPWIND_AVERAGE == YES
      alphaR = MAX(emf->SzR[k][j][i], emf->SzR[k][j][i+1]);
      alphaR = MAX(0,alphaR);
      alphaL = MIN(emf->SzL[k][j][i], emf->SzL[k][j][i+1]);
      alphaL = MIN(0,alphaL);
      
      scrh = 1.0/(alphaR - alphaL);
      wS =  alphaR*scrh;
      wN = -alphaL*scrh;
      aL = wS*emf->axL[k][j][i] + wN*emf->axL[k+1][j][i];
      aR = wS*emf->axR[k][j][i] + wN*emf->axR[k+1][j][i];
      dL = wS*emf->dxL[k][j][i] + wN*emf->dxL[k+1][j][i];
      dR = wS*emf->dxR[k][j][i] + wN*emf->dxR[k+1][j][i];
      #else
      aL = 0.5*(emf->axL[k][j][i] + emf->axL[k+1][j][i]);
      aR = 0.5*(emf->axR[k][j][i] + emf->axR[k+1][j][i]);
      dL = 0.5*(emf->dxL[k][j][i] + emf->dxL[k+1][j][i]);
      dR = 0.5*(emf->dxR[k][j][i] + emf->dxR[k+1][j][i]);
      #endif
      #endif

      phi = dR*bR[i] - dL*bL[i];
      emf->Ex2e[k][j][i] = (aL*vL[i]*bL[i] + aR*vR[i]*bR[i]) - phi;
    }
    #endif  /* INCLUDE_KDIR */
  }}
  #endif

/* --------------------------------------------------------
   X2. Reconstruct along the x2-direction:

     o -vy(k+1/2) and Bz(k+1/2) --> (i, j+1/2, k+1/2)
     o +vy(i+1/2) and Bx(i+1/2) --> (i+1/2, j+1/2, k)
   -------------------------------------------------------- */

  #if INCLUDE_JDIR
  for (k = emf->kbeg; k <= emf->kend; k++){ 
  for (i = emf->ibeg; i <= emf->iend; i++){
  
    #if INCLUDE_KDIR
    #if SHOCK_FLATTENING == MULTID
    for (j = emf->jbeg; j <= emf->jend+1; j++){
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k+1][j][i] & FLAG_MINMOD) ) flag[j] |= FLAG_MINMOD;
      else flag[j] = 0;
    }
    #endif
    ArrayReconstruct (emf->exk, flag, i, j, k, JDIR, vL, vR, recV, grid);  /* - vy */
    ArrayReconstruct (Bx3s,     flag, i, j, k, JDIR, bL, bR, recB, grid); 
    for (j = emf->jbeg; j <= emf->jend; j++){
      #if BACKGROUND_FIELD == YES
      BackgroundField (xr[i], yr[j], z[k], B0);
      bR[j] += B0[KDIR];
      bL[j] += B0[KDIR];
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      aL  = aR = 0.5;
      dL = 0.5*(emf->dyL[k][j][i] + emf->dyL[k+1][j][i]);
      dR = 0.5*(emf->dyR[k][j][i] + emf->dyR[k+1][j][i]);
      #else
      #if UPWIND_AVERAGE == YES
      alphaR = MAX(emf->SzR[k][j][i], emf->SzR[k][j+1][i]);
      alphaR = MAX(0,alphaR);
      alphaL = MIN(emf->SzL[k][j][i], emf->SzL[k][j+1][i]);
      alphaL = MIN(0,alphaL);
      
      scrh = 1.0/(alphaR - alphaL);
      wS =  alphaR*scrh;
      wN = -alphaL*scrh;
      aL = wS*emf->ayL[k][j][i] + wN*emf->ayL[k+1][j][i];
      aR = wS*emf->ayR[k][j][i] + wN*emf->ayR[k+1][j][i];
      dL = wS*emf->dyL[k][j][i] + wN*emf->dyL[k+1][j][i];
      dR = wS*emf->dyR[k][j][i] + wN*emf->dyR[k+1][j][i];
      #else
      aL = 0.5*(emf->ayL[k][j][i] + emf->ayL[k+1][j][i]);
      aR = 0.5*(emf->ayR[k][j][i] + emf->ayR[k+1][j][i]);
      dL = 0.5*(emf->dyL[k][j][i] + emf->dyL[k+1][j][i]);
      dR = 0.5*(emf->dyR[k][j][i] + emf->dyR[k+1][j][i]);
      #endif
      #endif

      phi = dR*bR[j] - dL*bL[j];
      emf->Ex1e[k][j][i] = (aL*vL[j]*bL[j] + aR*vR[j]*bR[j]) + phi;
    }
    #endif /* INCLUDE_KDIR */
    
    #if INCLUDE_IDIR
    #if SHOCK_FLATTENING == MULTID
    for (j = emf->jbeg; j <= emf->jend+1; j++){
      if (   (d->flag[k][j][i+1]   & FLAG_MINMOD)
          || (d->flag[k][j][i] & FLAG_MINMOD) ) flag[j] |= FLAG_MINMOD;
      else flag[j] = 0;
    }
    #endif
    ArrayReconstruct (emf->ezi, flag, i, j, k, JDIR, vL, vR, recV, grid);  /* + vy */
    ArrayReconstruct (Bx1s,     flag, i, j, k, JDIR, bL, bR, recB, grid);
    for (j = emf->jbeg; j <= emf->jend; j++){

      #if BACKGROUND_FIELD == YES
      BackgroundField (xr[i], yr[j], z[k], B0);
      bR[j] += B0[IDIR];
      bL[j] += B0[IDIR];
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      aL  = aR = 0.5;
      dL = 0.5*(emf->dyL[k][j][i] + emf->dyL[k][j][i+1]);
      dR = 0.5*(emf->dyR[k][j][i] + emf->dyR[k][j][i+1]);
      #else
      #if UPWIND_AVERAGE == YES
      alphaR = MAX(emf->SxR[k][j][i], emf->SxR[k][j+1][i]);
      alphaR = MAX(0,alphaR);
      alphaL = MIN(emf->SxL[k][j][i], emf->SxL[k][j+1][i]);
      alphaL = MIN(0,alphaL);
      
      scrh = 1.0/(alphaR - alphaL);
      wS =  alphaR*scrh;
      wN = -alphaL*scrh;
      aL = wS*emf->ayL[k][j][i] + wN*emf->ayL[k][j][i+1];
      aR = wS*emf->ayR[k][j][i] + wN*emf->ayR[k][j][i+1];
      dL = wS*emf->dyL[k][j][i] + wN*emf->dyL[k][j][i+1];
      dR = wS*emf->dyR[k][j][i] + wN*emf->dyR[k][j][i+1];
      #else
      aL = 0.5*(emf->ayL[k][j][i] + emf->ayL[k][j][i+1]);
      aR = 0.5*(emf->ayR[k][j][i] + emf->ayR[k][j][i+1]);
      dL = 0.5*(emf->dyL[k][j][i] + emf->dyL[k][j][i+1]);
      dR = 0.5*(emf->dyR[k][j][i] + emf->dyR[k][j][i+1]);
      #endif
      #endif

      phi = dR*bR[j] - dL*bL[j];
      emf->Ex3e[k][j][i] += (aL*vL[j]*bL[j] + aR*vR[j]*bR[j]) - phi;
    }
    #endif /* INCLUDE_IDIR */

  }}
  #endif

/* --------------------------------------------------------
   X3. Reconstruct along the x3-direction:

     o -vz(i+1/2) and Bx(i+1/2) --> (i+1/2, j, k+1/2)
     o +vz(j+1/2) and By(j+1/2) --> (i, j+1/2, k+1/2)
   -------------------------------------------------------- */

  #if INCLUDE_KDIR
  for (j = emf->jbeg; j <= emf->jend; j++){ 
  for (i = emf->ibeg; i <= emf->iend; i++){

    #if INCLUDE_IDIR
    #if SHOCK_FLATTENING == MULTID
    for (k = emf->kbeg; k <= emf->kend+1; k++){ 
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k][j][i+1] & FLAG_MINMOD) ) flag[k] |= FLAG_MINMOD;
      else flag[k] = 0;
    }
    #endif
    ArrayReconstruct (emf->eyi, flag, i, j, k, KDIR, vL, vR, recV, grid);  /* - vz */
    ArrayReconstruct (Bx1s,     flag, i, j, k, KDIR, bL, bR, recB, grid);  
    for (k = emf->kbeg; k <= emf->kend; k++){ 
      #if BACKGROUND_FIELD == YES
      BackgroundField (xr[i], y[j], zr[k], B0);
      bR[k] += B0[IDIR];
      bL[k] += B0[IDIR];
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      aL  = aR = 0.5;
      dL = 0.5*(emf->dzL[k][j][i] + emf->dzL[k][j][i+1]);
      dR = 0.5*(emf->dzR[k][j][i] + emf->dzR[k][j][i+1]);
      #else 
      #if UPWIND_AVERAGE == YES
      alphaR = MAX(emf->SxR[k][j][i], emf->SxR[k+1][j][i]);
      alphaR = MAX(0,alphaR);
      alphaL = MIN(emf->SxL[k][j][i], emf->SxL[k+1][j][i]);
      alphaL = MIN(0,alphaL);
      
      scrh = 1.0/(alphaR - alphaL);
      wS =  alphaR*scrh;
      wN = -alphaL*scrh;
      aL = wS*emf->azL[k][j][i] + wN*emf->azL[k][j][i+1];
      aR = wS*emf->azR[k][j][i] + wN*emf->azR[k][j][i+1];
      dL = wS*emf->dzL[k][j][i] + wN*emf->dzL[k][j][i+1];
      dR = wS*emf->dzR[k][j][i] + wN*emf->dzR[k][j][i+1];
      #else
      aL = 0.5*(emf->azL[k][j][i] + emf->azL[k][j][i+1]);
      aR = 0.5*(emf->azR[k][j][i] + emf->azR[k][j][i+1]);
      dL = 0.5*(emf->dzL[k][j][i] + emf->dzL[k][j][i+1]);
      dR = 0.5*(emf->dzR[k][j][i] + emf->dzR[k][j][i+1]);
      #endif
      #endif

      phi = dR*bR[k] - dL*bL[k];
      emf->Ex2e[k][j][i] += (aL*vL[k]*bL[k] + aR*vR[k]*bR[k]) + phi;
    }
    #endif /* INCLUDE_IDIR */

    #if INCLUDE_JDIR
    #if SHOCK_FLATTENING == MULTID
    for (k = emf->kbeg; k <= emf->kend+1; k++){ 
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k][j+1][i] & FLAG_MINMOD) ) flag[k] |= FLAG_MINMOD;
      else flag[k] = 0;
    }
    #endif
    ArrayReconstruct (emf->exj, flag, i, j, k, KDIR, vL, vR, recV, grid);   /* + vz */
    ArrayReconstruct (Bx2s,     flag, i, j, k, KDIR, bL, bR, recB, grid);  
    for (k = emf->kbeg; k <= emf->kend; k++){ 
      #if BACKGROUND_FIELD == YES
      BackgroundField (x[i], yr[j], zr[k], B0);
      bR[k] += B0[JDIR];
      bL[k] += B0[JDIR];
      #endif

      #if CT_EMF_AVERAGE == UCT_GFORCE
      aL  = aR = 0.5;
      dL = 0.5*(emf->dzL[k][j][i] + emf->dzL[k][j+1][i]);
      dR = 0.5*(emf->dzR[k][j][i] + emf->dzR[k][j+1][i]);
      #else
      #if UPWIND_AVERAGE == YES
      alphaR = MAX(emf->SyR[k][j][i], emf->SyR[k+1][j][i]);
      alphaR = MAX(0,alphaR);
      alphaL = MIN(emf->SyL[k][j][i], emf->SyL[k+1][j][i]);
      alphaL = MIN(0,alphaL);
      
      scrh = 1.0/(alphaR - alphaL);
      wS =  alphaR*scrh;
      wN = -alphaL*scrh;
      aL = wS*emf->azL[k][j][i] + wN*emf->azL[k][j+1][i];
      aR = wS*emf->azR[k][j][i] + wN*emf->azR[k][j+1][i];
      dL = wS*emf->dzL[k][j][i] + wN*emf->dzL[k][j+1][i];
      dR = wS*emf->dzR[k][j][i] + wN*emf->dzR[k][j+1][i];
      #else
      aL = 0.5*(emf->azL[k][j][i] + emf->azL[k][j+1][i]);
      aR = 0.5*(emf->azR[k][j][i] + emf->azR[k][j+1][i]);
      dL = 0.5*(emf->dzL[k][j][i] + emf->dzL[k][j+1][i]);
      dR = 0.5*(emf->dzR[k][j][i] + emf->dzR[k][j+1][i]);
      #endif
      #endif

      phi = dR*bR[k] - dL*bL[k];
      emf->Ex1e[k][j][i] += (aL*vL[k]*bL[k] + aR*vR[k]*bR[k]) - phi;
    }
    #endif /* INCLUDE_JDIR */
  }}
  #endif
}
#undef UPWIND_AVERAGE

/* ********************************************************************* */
void CT_EMF_Flux(const Data *d, const EMF *emf, Grid *grid)
/*!
 * Reconstruct electric field at faces to the edges.
 * The electric field interpolated from the 1D Riemann solver and
 * defined as E = E(pnt) + 2*E(diff), that is, the sum of the smooth
 * flux term plust \e twice the diffusive part.
 *********************************************************************** */
{
  int i, j, k;
#if RECONSTRUCTION == LINEAR
  #if LIMITER == DEFAULT
  int recE = MC_LIM;
  int recD = MC_LIM;
  #else
  int recE = LIMITER;
  int recD = LIMITER;
  #endif
#else
  int recE = RECONSTRUCTION;
  int recD = RECONSTRUCTION;
#endif
  DIM_EXPAND(double ***Bx1s = d->Vs[BX1s];  ,
             double ***Bx2s = d->Vs[BX2s];  ,
             double ***Bx3s = d->Vs[BX3s];)

  static char *flag;
  static double *eL, *eR, *dL, *dR;

  if (eL == NULL) {
    eL = ARRAY_1D(NMAX_POINT, double);
    eR = ARRAY_1D(NMAX_POINT, double);
    dL = ARRAY_1D(NMAX_POINT, double);
    dR = ARRAY_1D(NMAX_POINT, double);
    flag = ARRAY_1D(NMAX_POINT, char);
  }

/* --------------------------------------------------------
   X1. Reconstruct along the x1-direction:

     o Ezj(j+1/2) and Dy(j+1/2) --> (i+1/2, j+1/2, k)
     o Eyk(k+1/2) and Dz(k+1/2) --> (i+1/2, j, k+1/2)
   -------------------------------------------------------- */

  #if INCLUDE_IDIR
  for (k = emf->kbeg; k <= emf->kend; k++){ 
  for (j = emf->jbeg; j <= emf->jend; j++){

    #if INCLUDE_JDIR
    #if SHOCK_FLATTENING == MULTID
    for (i = emf->ibeg; i <= emf->iend+1; i++){
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k][j+1][i] & FLAG_MINMOD) ) flag[i] |= FLAG_MINMOD;
      else flag[i] = 0;
    }
    #endif
    ArrayReconstruct (emf->ezj,     flag, i, j, k, IDIR, eL, eR, recE, grid);
    ArrayReconstruct (emf->ezj_dff, flag, i, j, k, IDIR, dL, dR, recD, grid);    
    for (i = emf->ibeg; i <= emf->iend; i++){
      emf->Ex3e[k][j][i] = 0.25*(eL[i] + eR[i]) + 0.5*(dL[i] + dR[i]);
    }
    #endif  /* INCLUDE_JDIR */

    #if INCLUDE_KDIR
    #if SHOCK_FLATTENING == MULTID
    for (i = emf->ibeg; i <= emf->iend+1; i++){
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k+1][j][i] & FLAG_MINMOD) ) flag[i] |= FLAG_MINMOD;
      else flag[i] = 0;
    }
    #endif
    ArrayReconstruct (emf->eyk,     flag, i, j, k, IDIR, eL, eR, recE, grid);
    ArrayReconstruct (emf->eyk_dff, flag, i, j, k, IDIR, dL, dR, recD, grid);
    for (i = emf->ibeg; i <= emf->iend; i++){
      emf->Ex2e[k][j][i] = 0.25*(eL[i] + eR[i]) + 0.5*(dL[i] + dR[i]);
    }
    #endif  /* INCLUDE_KDIR */

  }}
  #endif

/* --------------------------------------------------------
   X2. Reconstruct along the x2-direction:

     o Exk(k+1/2) and Dz(k+1/2) --> (i, j+1/2, k+1/2)
     o Ezi(i+1/2) and Dx(i+1/2) --> (i+1/2, j+1/2, k)
   -------------------------------------------------------- */

  #if INCLUDE_JDIR
  for (k = emf->kbeg; k <= emf->kend; k++){ 
  for (i = emf->ibeg; i <= emf->iend; i++){

    #if INCLUDE_KDIR
    #if SHOCK_FLATTENING == MULTID
    for (j = emf->jbeg; j <= emf->jend+1; j++){
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k+1][j][i] & FLAG_MINMOD) ) flag[j] |= FLAG_MINMOD;
      else flag[j] = 0;
    }
    #endif
    ArrayReconstruct (emf->exk,     flag, i, j, k, JDIR, eL, eR, recE, grid);
    ArrayReconstruct (emf->exk_dff, flag, i, j, k, JDIR, dL, dR, recD, grid);    
    for (j = emf->jbeg; j <= emf->jend; j++){
      emf->Ex1e[k][j][i] = 0.25*(eL[j] + eR[j]) + 0.5*(dL[j] + dR[j]);
    }
    #endif  /* INCLUDE_KDIR */

    #if INCLUDE_IDIR
    #if SHOCK_FLATTENING == MULTID
    for (j = emf->jbeg; j <= emf->jend+1; j++){
      if (   (d->flag[k][j][i+1]   & FLAG_MINMOD)
          || (d->flag[k][j][i] & FLAG_MINMOD) ) flag[j] |= FLAG_MINMOD;
      else flag[j] = 0;
    }
    #endif
    ArrayReconstruct (emf->ezi,     flag, i, j, k, JDIR, eL, eR, recE, grid);
    ArrayReconstruct (emf->ezi_dff, flag, i, j, k, JDIR, dL, dR, recD, grid);
    for (j = emf->jbeg; j <= emf->jend; j++){
      emf->Ex3e[k][j][i] += 0.25*(eL[j] + eR[j]) + 0.5*(dL[j] + dR[j]);
    }
    #endif /* INCLUDE_IDIR */

  }}
  #endif

/* --------------------------------------------------------
   X3. Reconstruct along the x3-direction:

     o Eyi(i+1/2) and Dx(i+1/2) --> (i+1/2, j, k+1/2)
     o Exj(j+1/2) and Dy(j+1/2) --> (i, j+1/2, k+1/2)
   -------------------------------------------------------- */

  #if INCLUDE_KDIR
  for (j = emf->jbeg; j <= emf->jend; j++){ 
  for (i = emf->ibeg; i <= emf->iend; i++){

    #if INCLUDE_IDIR
    #if SHOCK_FLATTENING == MULTID
    for (k = emf->kbeg; k <= emf->kend+1; k++){ 
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k][j][i+1] & FLAG_MINMOD) ) flag[k] |= FLAG_MINMOD;
      else flag[k] = 0;
    }
    #endif
    ArrayReconstruct (emf->eyi,     flag, i, j, k, KDIR, eL, eR, recE, grid);
    ArrayReconstruct (emf->eyi_dff, flag, i, j, k, KDIR, dL, dR, recD, grid);
    for (k = emf->kbeg; k <= emf->kend; k++){
      emf->Ex2e[k][j][i] += 0.25*(eL[k] + eR[k]) + 0.5*(dL[k] + dR[k]);
    }
    #endif /* INCLUDE_IDIR */

    #if INCLUDE_JDIR
    #if SHOCK_FLATTENING == MULTID
    for (k = emf->kbeg; k <= emf->kend+1; k++){ 
      if (   (d->flag[k][j][i]   & FLAG_MINMOD)
          || (d->flag[k][j+1][i] & FLAG_MINMOD) ) flag[k] |= FLAG_MINMOD;
      else flag[k] = 0;
    }
    #endif
    ArrayReconstruct (emf->exj,     flag, i, j, k, KDIR, eL, eR, recE, grid);
    ArrayReconstruct (emf->exj_dff, flag, i, j, k, KDIR, dL, dR, recD, grid);
    for (k = emf->kbeg; k <= emf->kend; k++){
      emf->Ex1e[k][j][i] += 0.25*(eL[k] + eR[k]) + 0.5*(dL[k] + dR[k]);
    }
    #endif /* INCLUDE_IDIR */

  }}
  #endif
 
}
