/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Functions for computing/retrieving the mean aziumthal velocity.

  \authors A. Mignone (mignone@to.infn.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)
         
  \date    Feb 10, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/*! Defines a 2D array containing the azimuthally-averaged velocity.
    When the average is done along the X2 direction (Cartesian and polar
    geometries) the array should have size (NX3_TOT, NX1_TOT).
    When the average is performed along the X3 direction (spherical geometry)
    the array has dimensions (NX2_TOT, NX1_TOT)                           */
static double **wA;  
/*static double **wAx1, **wAx2, **wAx3;*/  /* -- average orbital speed at
                                            x1, x2 and x3 faces -- */
                                            
/* ********************************************************************* */
void FARGO_Initialize(void)
/*!
 * Initialize FARGO module.
 *
 *********************************************************************** */
{
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
  if (wA == NULL) wA = ARRAY_2D(NX3_TOT, NX1_TOT, double);
  #else
  if (wA == NULL) wA = ARRAY_2D(NX2_TOT, NX1_TOT, double);
  #endif
}

/* ********************************************************************* */
void FARGO_AverageVelocity(const Data *d, Grid *grid)
/*!
 * Compute the background orbital velocity as a 2D array by averaging
 * the total azimuthal velocity along the orbital coordinate y/phi/phi.
 * Update the residual velocity as well.
 *
 *********************************************************************** */
{
#if FARGO_NSTEP_AVERAGE > 0
  int i,j,k;
  double th, w;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  #if GEOMETRY == CARTESIAN
  double ***vphi = d->Vc[VX2];
  #else
  double ***vphi = d->Vc[iVPHI];
  #endif

  if ( (g_stepNumber%FARGO_NSTEP_AVERAGE == 0) && g_stepNumber != 0){

    #ifdef PARALLEL
    double w_recv=0.0, w_sum=0.0;
    int count, coords[3], src, dst;
    MPI_Comm cartcomm;
    MPI_Status status;
    #endif
     
  /* ------------------------------------------------------
     1. Fill ghost zones
     ------------------------------------------------------ */

    Boundary(d, ALL_DIR, grid); 

  /* ------------------------------------------------------
     2. get ranks of the upper and lower procs
     ------------------------------------------------------ */

    #ifdef PARALLEL
    AL_Get_cart_comm(SZ, &cartcomm);
    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[SDIR] += 1;
    MPI_Cart_rank(cartcomm, coords, &dst);

    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[SDIR] -= 1;
    MPI_Cart_rank(cartcomm, coords, &src);
    #endif

  /* ------------------------------------------------------
     3. Compute average orbital speed
        in CART. or POLAR geometries
     ------------------------------------------------------ */

    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    KTOT_LOOP(k) ITOT_LOOP(i) {

    /* ------------------------------------------
       3a. Compute average orbital speed
       ------------------------------------------ */

      w = 0.0;
      JDOM_LOOP(j) w += wA[k][i] + vphi[k][j][i]; /* Sum total velocity */
      
      if (grid->nproc[JDIR] > 1){
        #ifdef PARALLEL
        w_sum = w;
        for (count=1; count < grid->nproc[JDIR]; count++ ){
          MPI_Sendrecv(&w, 1, MPI_DOUBLE, dst, 0,  &w_recv, 1,
                              MPI_DOUBLE, src, 0, cartcomm, &status);
          w      = w_recv;
          w_sum += w;
        }
        w = w_sum/(double)(grid->np_int_glob[JDIR]);
        #endif
      }else{
         w = w/(double)NX2;
      }

    /* ------------------------------------------
       3b. Update residual
       ------------------------------------------ */

      JDOM_LOOP(j) vphi[k][j][i] += wA[k][i] - w;
      wA[k][i] = w;
    }
    #endif

  /* ------------------------------------------------------
     4. Compute average orbital speed
        in SPHERICAL geom.
     ------------------------------------------------------ */

    #if GEOMETRY == SPHERICAL
    JTOT_LOOP(j) ITOT_LOOP(i) {

    /* ------------------------------------------
       4a. Compute average orbital speed
       ------------------------------------------ */

      w = 0.0;
      KDOM_LOOP(k) w += wA[j][i] + vphi[k][j][i];
      if (grid->nproc[KDIR] > 1){
        #ifdef PARALLEL
        w_sum = w;
        for (count=1; count < grid->nproc[KDIR]; count++ ){
          MPI_Sendrecv(&w, 1, MPI_DOUBLE, dst, 0,  &w_recv, 1,
                              MPI_DOUBLE, src, 0, cartcomm, &status);
          w      = w_recv;
          w_sum += w;
        }
        w = w_sum/(double)(grid->np_int_glob[KDIR]);
        #endif
      }else{
         w = w/(double)NX3;
      }

    /* ------------------------------------------
       4b. Update residual
       ------------------------------------------ */
    
      KDOM_LOOP(k) vphi[k][j][i] += wA[j][i] - w;
      wA[j][i] = w;
    }
    #endif
  } 
#endif /* FARGO_NSTEP_AVERAGE > 0 */
}

/* ********************************************************************* */
void FARGO_ComputeTotalVelocity(const Data *d, double ***vtot, Grid *grid)
/*!
 * Add the mean background contribution to the residual
 * velocity in order to obtain the total velocity.
 *
 * \param [in,out] d      pointer to PLUTO Data structure
 * \param [in,out] vtot   array where the total velocity is stored.
 *                        (can be also d->Vc[VX2])
 * \param [in]     grid   pointer to array of Grid structures
 *
 * \return This function has no return value.
 *********************************************************************** */
{
  int i,j,k;
#if GEOMETRY == CARTESIAN
  double ***vphi = d->Vc[VX2];
#else
  double ***vphi = d->Vc[iVPHI];
#endif

/* --------------------------------------------------------
   2a. Add background shear velocity to the fluid
   -------------------------------------------------------- */

  TOT_LOOP(k,j,i){
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    vtot[k][j][i] = wA[k][i] + vphi[k][j][i];
    #elif GEOMETRY == SPHERICAL
    vtot[k][j][i] = wA[j][i] + vphi[k][j][i];
    #else
    printLog ("! FARGO not supported in this geometry\n");
    QUIT_PLUTO(1);
    #endif 
  }

/* --------------------------------------------------------
   2a. Add background shear velocity to particles.
       Note: we do not call Particles_Interpolate since
       the array to be interpolated has only 2 dimensions
       and not three.
   -------------------------------------------------------- */

#if PARTICLES
  int i1,j1,k1;
  particleNode *curNode;
  Particle *p;
  double v;
  static double ***wp;

  if (wp == NULL)  wp = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);

  PARTICLES_LOOP(curNode, d->PHead){
    p = &(curNode->p);

  /* -- Interpolate fargo shift velocity at particle position -- */

    Particles_GetWeights(p, p->cell, wp, grid);  
    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];

    v = 0.0;
    for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
    for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
    for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
      v += wp[k1][j1][i1]*wA[k + k1][i + i1];
    }}}
/*  v = -SB_Q*SB_OMEGA*p->coord[IDIR]; */
    p->speed[JDIR] += v;
  }
#endif

}

/* ********************************************************************* */
double **FARGO_Velocity(void)
/*!
 * Return a pointer to the background orbital velocity ::wA.
 *
 *********************************************************************** */
{
  return wA;
}
/* ********************************************************************* */
void FARGO_Restart(const Data *d, Grid *grid)
/*!
 * Re-compute residual and average orbital speed after restarting from file.
 * Warning: this operation is not possible when the background velocity is
 *          computed dynamically and the output contains the residual.
 *          The following table summarizes the possible combination.
 *

  <CENTER>
  FARGO_NSTEP_AVERAGE | FARGO_OUTPUT_VTOT |  Restart_Possible ? 
  --------------------|-------------------|---------------------------
       >0             |   NO  (v')        |  No
       >0             |   YES (v)         |  Yes (w = <v>, v' = v-w)
       -1             |   NO  (v')        |  Yes  
       -1             |   YES (v)         |  Yes (v' = v - w)
  </CENTER>

 *********************************************************************** */
{
  int i,j,k;

#if FARGO_NSTEP_AVERAGE > 0
  
  #if FARGO_OUTPUT_VTOT == NO
  printLog ("! FARGO_Restart(): cannot compute background velocity since\n");
  printLog ("                   data contains the residual only.\n");
  QUIT_PLUTO(1);
  #endif
  
  #if (GEOMETRY == CARTESIAN) || (GEOMETRY == POLAR) 
  KTOT_LOOP(k) ITOT_LOOP(i) wA[k][i] = 0.0;
  #elif GEOMETRY == SPHERICAL
  JTOT_LOOP(j) ITOT_LOOP(i) wA[j][i] = 0.0;
  #endif
  
  FARGO_AverageVelocity(d, grid);

#else

  #if FARGO_OUTPUT_VTOT == YES
  #if (GEOMETRY == CARTESIAN) || (GEOMETRY == POLAR) 
  TOT_LOOP(k,j,i) d->Vc[VX2][k][j][i] -= wA[k][i]; 
  #elif GEOMETRY == SPHERICAL
  TOT_LOOP(k,j,i) d->Vc[VX3][k][j][i] -= wA[j][i];
  #endif
  #endif
  
#endif

}
