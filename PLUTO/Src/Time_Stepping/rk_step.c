/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advance equations with Runge Kutta time integrators.

  Main driver for RK split/unsplit integrations and finite difference
  methods (RK3).
  Time stepping include Euler, RK2 and RK3.

  \authors A. Mignone (mignone@to.infn.it)\n
  \date    Nov 11, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* Weight factor for 2nd stage of RK integrators */

#if TIME_STEPPING == RK2
 #define w0  0.5
 #define wc  0.5
#elif TIME_STEPPING == RK3
 #define w0 0.75
 #define wc 0.25
#endif

static void SolutionCorrect(Data *, timeStep *, Data_Arr, Data_Arr, double, Grid *);

/* ********************************************************************* */
int AdvanceStep (Data *d, timeStep *Dts, Grid *grid)
/*!
 * Advance the equations by a single time step using unsplit 
 * integrators based on the method of lines.
 *
 * \param [in,out]      d  pointer to Data structure
 * \param [in,out]    Dts  pointer to time step structure
 * \param [in]       grid  pointer to array of Grid structures
 *    
 *********************************************************************** */
{
  int  i, j, k, nv;
  static double  one_third = 1.0/3.0;
  static Data_Arr U0;
  static double ***Bs0[3];
  RBox   box;
#if PARTICLES
  Data_Arr Vpnt;
  static Data_Arr Vhalf;
#endif

  RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &box);

  #if RADIATION && RADIATION_IMEX_SSP2 && (TIME_STEPPING != RK2)
  print ("! AdvanceStep(): RADIATION_IMEX_SSP2 requires TIME_STEPPING == RK2\n");
  QUIT_PLUTO(1);
  #endif

  #if RADIATION && RADIATION_IMEX_SSP2
  static Data_Arr Srad1, Srad2;
  static double rk_gamma = 0.29289321881;
  double dt_rk1, dt_rk2 ,dt_rk3 ;

  dt_rk1 = g_dt*rk_gamma ;
  dt_rk2 = g_dt*(1.0-3.0*rk_gamma);
  dt_rk3 = g_dt*0.5*(1.0-rk_gamma);
  #endif

/* --------------------------------------------------------
   0. Allocate memory 
   -------------------------------------------------------- */

  if (U0 == NULL){
    U0 = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    #if PARTICLES
    Vhalf = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif

    #ifdef STAGGERED_MHD
    DIM_EXPAND(
      Bs0[IDIR] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
      Bs0[JDIR] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
      Bs0[KDIR] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    )
    #endif

    #if RADIATION
    #if RADIATION_IMEX_SSP2
    Srad1 = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, 4, double);
    Srad2 = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, 4, double);
    #endif
    #if RADIATION_VAR_OPACITIES == NO
    g_totalOpacity = g_absorptionCoeff + g_scatteringCoeff ;
    #endif
    #endif
  }

/* --------------------------------------------------------
   1. Predictor step (EULER, RK2, RK3)

      After baoundaries have been set we flag zones lying 
      in a shock. 
      This is useful for shock flattening or 
      entropy/energy selective update.
   -------------------------------------------------------- */

/* -- 1a. Set boundary conditions -- */

  g_intStage = 1;
  #if RING_AVERAGE > 1
  PrimToCons3D (d->Vc, d->Uc, &box);
  RingAverageCons(d, grid);
  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);
  #endif
  Boundary (d, ALL_DIR, grid);
  #if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH) 
  FlagShock (d, grid);
  #endif

/* -- 1b. Convert primitive to conservative, save initial stage  -- */

  PrimToCons3D(d->Vc, d->Uc, &box);
  RBoxCopy (&box, U0, d->Uc, NVAR, CONS_ARRAY);
#ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i) Bs0[nv][k][j][i] = d->Vs[nv][k][j][i];
#endif

/* -- 1c. Compute Particles feedback / Advance with pred. step -- */

  #if PARTICLES == PARTICLES_CR
  Particles_CR_ComputeForce (d->Vc, d, grid);
  NVAR_LOOP(nv) TOT_LOOP(k,j,i) Vhalf[nv][k][j][i] = 0.5*d->Vc[nv][k][j][i];
  #elif PARTICLES == PARTICLES_DUST
  Particles_Dust_ComputeForce (d->Vc, d, grid);
  NVAR_LOOP(nv) TOT_LOOP(k,j,i) Vhalf[nv][k][j][i] = 0.5*d->Vc[nv][k][j][i];
  //#elif PARTICLES == PARTICLES_LP
  //Particles_LP_Update (d, Dts, g_dt, grid);
  #endif

/* -- 1d. Advance conservative variables array -- */

#if RADIATION
  #if RADIATION_IMEX_SSP2
  RadStep3D (d->Uc, d->Vc, Srad1, d->flag, &box, dt_rk1);
  Boundary (d, ALL_DIR, grid); 
  #endif
#endif

/* CheckData (d, grid, "Before Predictor"); */
  UpdateStage(d, d->Uc, d->Vs, NULL, g_dt, Dts, grid);
  #if TIME_STEP_CONTROL == YES
  SolutionCorrect(d, Dts, U0, Bs0, g_dt, grid);
  #endif

  #if RING_AVERAGE > 1
  RingAverageCons(d, grid);
  #endif

#if RADIATION
  #if RADIATION_IMEX_SSP2
  AddRadSource1(d->Uc, Srad1, &box, dt_rk2);
  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);
  RadStep3D (d->Uc, d->Vc, Srad2, d->flag, &box, dt_rk1);
  #else
  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);
  RadStep3D (d->Uc, d->Vc, NULL, d->flag, &box, g_dt);
  #endif 
#endif

#ifdef STAGGERED_MHD
  CT_AverageStaggeredFields (d->Vs, d->Uc, &box, grid);
#endif
  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);

/* -- 1e. Advance particles (w/ feedback) by a full step -- */

#if   (PARTICLES == PARTICLES_CR   && PARTICLES_CR_FEEDBACK == YES)   \
   || (PARTICLES == PARTICLES_DUST && PARTICLES_DUST_FEEDBACK == YES) 
  NVAR_LOOP(nv) TOT_LOOP(k,j,i) Vhalf[nv][k][j][i] += 0.5*d->Vc[nv][k][j][i];
  Vpnt  = d->Vc;  /* Save pointer */
  d->Vc = Vhalf;
  #if PARTICLES == PARTICLES_CR
  Particles_CR_Update(d, Dts, g_dt, grid);
  #elif PARTICLES == PARTICLES_DUST
  Particles_Dust_Update(d, Dts, g_dt, grid);
  #endif
  d->Vc = Vpnt;   /* Restore Pointer */
#endif  /* PARTICLES_XX_FEEDBACK */
  
/* --------------------------------------------------------
   2. Corrector step (RK2, RK3)
   -------------------------------------------------------- */

#if (TIME_STEPPING == RK2) || (TIME_STEPPING == RK3)

/* -- 2a. Set boundary conditions -- */

  g_intStage = 2;
  Boundary (d, ALL_DIR, grid);

/* -- 2b. Advance paticles & solution array -- */
  
  #if PARTICLES == PARTICLES_LP
  Particles_LP_Update (d, Dts, g_dt, grid);
  #endif
  UpdateStage(d, d->Uc, d->Vs, NULL, g_dt, Dts, grid);
  
  #if RADIATION && !RADIATION_IMEX_SSP2
  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);
  RadStep3D (d->Uc, d->Vc, NULL, d->flag, &box, g_dt);
  #endif 

  DOM_LOOP(k, j, i) NVAR_LOOP(nv){
    d->Uc[k][j][i][nv] = w0*U0[k][j][i][nv] + wc*d->Uc[k][j][i][nv];
  }
  #if RING_AVERAGE > 1
  RingAverageCons(d, grid);
  #endif

#if RADIATION
  #if RADIATION_IMEX_SSP2
  AddRadSource2(d->Uc, Srad1, Srad2, &box, dt_rk1, dt_rk3);
  #endif 
#endif

  #ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i) {
    d->Vs[nv][k][j][i] = w0*Bs0[nv][k][j][i] + wc*d->Vs[nv][k][j][i];
  }
  CT_AverageStaggeredFields (d->Vs, d->Uc, &box, grid);
  #endif

/* -- 2e. Apply FARGO orbital shift -- */

  #if (defined FARGO) && (TIME_STEPPING == RK2)
  FARGO_ShiftSolution (d->Uc, d->Vs, grid);
  #endif

/* -- 2f. Convert to Primitive -- */

  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);

#endif  /* TIME_STEPPING == RK2/RK3 */

/* --------------------------------------------------------
   3. Last corrector step (RK3 only) 
   -------------------------------------------------------- */

#if TIME_STEPPING == RK3
  #if (PARTICLES == PARTICLES_CR   && PARTICLES_CR_FEEDBACK == YES)   \
   || (PARTICLES == PARTICLES_DUST && PARTICLES_DUST_FEEDBACK == YES) 
  print ("! AdvanceStep(): RK3 algorithm not permitted with particles and feedback\n");
  QUIT_PLUTO(1);
  #endif

/* -- 3a. Set Boundary conditions -- */

  g_intStage = 3;
  Boundary (d, ALL_DIR, grid);

/* -- 3b. Update solution array -- */

  #if PARTICLES == PARTICLES_LP
  Particles_LP_Update (d, Dts, g_dt, grid);
  #endif
  UpdateStage(d, d->Uc, d->Vs, NULL, g_dt, Dts, grid);

  #if RADIATION
  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);
  RadStep3D (d->Uc, d->Vc, NULL, d->flag, &box, g_dt);
  #endif

  DOM_LOOP(k,j,i) NVAR_LOOP(nv){
    d->Uc[k][j][i][nv] = one_third*(U0[k][j][i][nv] + 2.0*d->Uc[k][j][i][nv]);
  }
  #if RING_AVERAGE > 1
  RingAverageCons(d, grid);
  #endif

  #ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i){
    d->Vs[nv][k][j][i] = (Bs0[nv][k][j][i] + 2.0*d->Vs[nv][k][j][i])/3.0;
  }
  CT_AverageStaggeredFields (d->Vs, d->Uc, &box, grid);
  #endif

/* -- 3c. Apply FARGO orbital shift -- */

  #ifdef FARGO
  FARGO_ShiftSolution (d->Uc, d->Vs, grid);
  #endif
  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);
#endif /* TIME_STEPPING == RK3 */

/* --------------------------------------------------------
   4. Particles update (no feedback)
   -------------------------------------------------------- */

  #if (PARTICLES == PARTICLES_CR   && PARTICLES_CR_FEEDBACK == NO)   \
   || (PARTICLES == PARTICLES_DUST && PARTICLES_DUST_FEEDBACK == NO) 
  NVAR_LOOP(nv) DOM_LOOP(k,j,i) Vhalf[nv][k][j][i] += 0.5*d->Vc[nv][k][j][i];
  Vpnt  = d->Vc;  /* Save pointer */
  d->Vc = Vhalf;
  #if PARTICLES == PARTICLES_CR
  Particles_CR_Update(d, Dts, g_dt, grid);
  #elif PARTICLES == PARTICLES_DUST
  Particles_Dust_Update(d, Dts, g_dt, grid);
  #endif
  d->Vc = Vpnt;   /* Restore Pointer */
#endif  /* PARTICLES */

/* --------------------------------------------------------
   5. Inject particles or update spectra
   -------------------------------------------------------- */
  
#if PARTICLES
  #if (PARTICLES == PARTICLES_LP) && (PARTICLES_LP_SPECTRA == YES)
  Boundary (d, ALL_DIR, grid); /* Spectra update requires interpolating
                                * fluid quantities at particle position. */
  Particles_LP_UpdateSpectra (d, g_dt, grid);
  #endif
  Particles_Inject(d,grid);
#endif

  return 0; /* -- step has been achieved, return success -- */
}

#if TIME_STEP_CONTROL == YES
/* ********************************************************************* */
void SolutionCorrect(Data *d, timeStep *Dts,
                     Data_Arr U0, Data_Arr Vs0, double dt0, Grid *grid)
/*
 * Recompute the time step (dt) based on most recent call and compare it
 * with the actual time step dt0 being used for the present update.
 * If dt < rmax*dt0 then lower the time step to dt = rsafe*dt.
 *
 * This function should be called after the first predictor step in
 * RK time-stepping function.
 *
 *********************************************************************** */
{
  int i,j,k,nv;
  Runtime *runtime = RuntimeGet();
  double dt = NextTimeStep(Dts, runtime, grid);
  double rmax  = 0.5;  // 0.65
  double rsafe = 1.0;  // 0.9 

  if (dt < rmax*dt0){
    dt = rsafe*dt;
    print ("! SolutionCorrect(): time step must be lowered (dt/dt0 = %f)\n", dt/dt0);
    DOM_LOOP(k,j,i){
      NVAR_LOOP(nv) {
        double dU = (d->Uc[k][j][i][nv] - U0[k][j][i][nv]);
        d->Uc[k][j][i][nv] = U0[k][j][i][nv] + dU*dt/dt0;
      }
    }
    DIM_LOOP(nv) {
      TOT_LOOP(k,j,i){
        double dV = d->Vs[nv][k][j][i] - Vs0[nv][k][j][i];
        d->Vs[nv][k][j][i] = Vs0[nv][k][j][i] + dV*dt/dt0;
      }
    }
  }

  g_dt = dt;
}
#endif
