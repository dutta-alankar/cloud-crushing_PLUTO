/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Integrate cooling and reaction source terms.

  Solve the system of ordinary differential equations describing
  optically thin cooling and ionization network (if any).
  We use an adaptive step size integration method and follow the 
  strategy outlined in section 5.1 of Tesileanu et al. (2008) for handling
  stiff cells.

  On output, a time-step estimate for the next time level is computed using
  the relative or absolute variation obtained during the integration of the ODE 
  (from t(n) --> t(n+1))  
  \f[
     \Delta t_c = \min_{ijk}\left[\frac{\Delta t^n M_R}{\epsilon}\right]
     \,\quad\rm{where}\quad
     \epsilon = \max\left(\left|\frac{p^{n+1}}{p^n} - 1\right|,\,
                |X^{n+1}-X^n|\right)
  \f]
  where \f$ M_R \f$ is the maximum cooling rate (defined by the global variable  
  ::g_maxCoolingRate) and X are the chemical species.
  
  \b References
     - "Simulating radiative astrophysical flows with the PLUTO code:
        a non-equilibrium, multi-species cooling function" \n
       Tesileanu, Mignone and Massaglia, A&A (2008) 488, 429

  \authors A. Mignone (mignone@to.infn.it)\n
           O. Tesileanu
           B. Vaidya
           A. Dutta (\date   Nov 24, 2021 : Added Townsend cooling)
  \date    April 09, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "cooling.h"

#if COOLING == GRACKLE
/* ********************************************************
    Global coordinates available in the Radiat() function
   ******************************************************** */

double gCooling_x1, gCooling_x2, gCooling_x3;

/* ********************************************************************* */
void CoolingSource (const Data *d, double dt, timeStep *Dts, Grid *grid)
/*!
 * Integrate cooling and reaction source terms.
 *
 * \param [in,out]  d     pointer to Data structure
 * \param [in]     dt     the time step to be taken
 * \param [out]    Dts    pointer to the Time_Step structure
 * \param [in]     grid   pointer to an array of Grid structures
 *
 *********************************************************************** */
{
    call_grackle(d, dt, Dts, grid, 0, 0, 0, 0);
}
#elif (COOLING != TOWNSEND) && (COOLING ! =GRACKLE)
/* ********************************************************
    Global coordinates available in the Radiat() function
   ******************************************************** */

double gCooling_x1, gCooling_x2, gCooling_x3;

/* ********************************************************************* */
void CoolingSource (const Data *d, double dt, timeStep *Dts, Grid *grid)
/*!
 * Integrate cooling and reaction source terms.
 *
 * \param [in,out]  d     pointer to Data structure
 * \param [in]     dt     the time step to be taken
 * \param [out]    Dts    pointer to the Time_Step structure
 * \param [in]     grid   pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int  nv, k, j, i, stiff, status;
  double err, scrh, min_tol = 2.e-5;
  double mu0, T0, T1, mu1, prs;
  double v0[NVAR_COOLING], v1[NVAR_COOLING], k1[NVAR_COOLING];
  double maxrate;
  intList var_list;
 
/* --------------------------------------------------------
   0. Create list of time-dependent variables for this
      step. Initialize coefficients.
   -------------------------------------------------------- */

  var_list.nvar    = NIONS+1;
  var_list.indx[0] = RHOE;
  for (nv = 0; nv < NIONS; nv++) var_list.indx[nv+1] = NFLX+nv;

  for (nv = NVAR_COOLING; nv--;  ) {
    v0[nv] = v1[nv] = k1[nv] = 0.0;
  }

/* -------------------------------------------------------- 
   1. Main loop on interior computational zones
   --------------------------------------------------------  */

  DOM_LOOP(k,j,i){
  
  /* ------------------------------------------------------
     1A. Define global coordinates.
         (useful for spatial-dependent cooling
     ------------------------------------------------------ */ 

    gCooling_x1 = grid->x[IDIR][i];
    gCooling_x2 = grid->x[JDIR][j];
    gCooling_x3 = grid->x[KDIR][k];

  /* ------------------------------------------------------
     1B. Skip integration if cell has been tagged with 
      FLAG_INTERNAL_BOUNDARY or FLAG_SPLIT_CELL 
      (AMR only).
     ------------------------------------------------------ */ 

    #if INTERNAL_BOUNDARY == YES
    if (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY) continue;
    #endif
    if (d->flag[k][j][i] & FLAG_SPLIT_CELL) continue;
    NVAR_LOOP(nv) v0[nv] = v1[nv] = d->Vc[nv][k][j][i];

  /* ------------------------------------------------------
     1C. Compute temperature and internal energy from
         density, pressure and concentrations.
     ------------------------------------------------------ */
    
    prs = v0[PRS];
    #if COOLING != TABULATED
    mu0 = MeanMolecularWeight(v0);
    #else
    double dummy0[4];
    mu0 = MeanMolecularWeight(v0, dummy0);
    #endif
    T0  = v0[PRS]/v0[RHO]*KELVIN*mu0;

    #if EOS == IDEAL
    v0[RHOE] = v1[RHOE] = prs/(g_gamma-1.0);
    #else
    v1[RHOE] = v0[RHOE] = InternalEnergy(v0, T0);
    #endif

    if (T0 <= 0.0){
      printLog ("! CoolingSource() line 123: negative initial temperature\n");
      printLog (" %12.6e  %12.6e\n",v0[RHOE], v0[RHO]);
      printLog (" at: %f %f\n",gCooling_x1, gCooling_x2);
      QUIT_PLUTO(1);
    }

  /* ------------------------------------------------------
     1D. Integrate the ODEs from t1 = 0...dt with
         increments dt1.

         - if the rhs is not stiff, try with a lower order
           method (e.g. RK2);
         - otherwise switch to Cash-Karp 5th order with
           adaptive step size dt1.
     ------------------------------------------------------ */

    Radiat(v0, k1);

    stiff = (fabs(k1[RHOE]) > 0.5/dt ? 1:0);
    if (!stiff){
      err  = SolveODE_RKF12 (v0, k1, v1, dt, &var_list);
      err /= min_tol;
    }

    if (stiff || err > 1.0){  /* Adaptive step size */
      double dt1, dts=dt;
      double t1   = 0.0;  /* t1 = 0, ..., dt */ 
      int    done = 0;

      dt1 = 0.5/(fabs(k1[RHOE]) + 1.0/dt);   /* Rough estimate to initial dt1 */
      do {
        
      /* -- Adjust time step dt1 so we do not overshoot dt -- */

        if ( (t1+dt1) >= dt) {
          dt1 = dt - t1;
          done = 1;
        }

        dts = SolveODE_CK45 (v0, k1, v1, dt1, min_tol, &var_list);
        t1 += dt1;  
        dt1 = dts;

        if (done) break;
        for (nv = NVAR_COOLING; nv--;  ) v0[nv] = v1[nv];
        Radiat(v0, k1);

      }while (!done);
    }

  /* ------------------------------------------------------
     1E. As asafety measure, constrain the ions to
         lie in [0,1] and pressure to be positive.
     ------------------------------------------------------ */

    NIONS_LOOP(nv){
      v1[nv] = MAX(v1[nv], 0.0);
      v1[nv] = MIN(v1[nv], 1.0);
    }
    #if COOLING == H2_COOL
     v1[X_H2] = MIN(v1[X_H2], 0.5);
    #endif
    
  /* -- pressure must be positive -- */
  
    #if COOLING != TABULATED
    mu1 = MeanMolecularWeight(v1);
    #else
    double dummy1[4];
    mu1 = MeanMolecularWeight(v1, dummy1);
    #endif
    
    #if EOS == IDEAL
    prs = v1[RHOE]*(g_gamma - 1.0);
    T1  = prs/v1[RHO]*KELVIN*mu1;
    #elif EOS == PVTE_LAW
    status = GetEV_Temperature(v1[RHOE], v1, &T1);
    prs    = v1[RHO]*T1/(KELVIN*mu1);
    #endif

    if (prs < 0.0) prs = g_smallPressure;

  /* ------------------------------------------------------
     1F. If gas cooled below lower threshold, redefine
         pressure so that T = g_minCoolingTemp.
     ------------------------------------------------------ */

    if (T1 < g_minCoolingTemp && T0 > g_minCoolingTemp)
      prs = g_minCoolingTemp*v1[RHO]/(KELVIN*mu1);

  /* ------------------------------------------------------
     1G. Suggest next time step based on 
         fractional variaton.
     ------------------------------------------------------ */

    err = fabs(prs/d->Vc[PRS][k][j][i] - 1.0);

    #if COOLING == MINEq
    for (nv = NFLX; nv < NFLX + NIONS - Fe_IONS; nv++) err = MAX(err, fabs(d->Vc[nv][k][j][i] - v1[nv]));
    #else
    NIONS_LOOP(nv) err = MAX(err, fabs(d->Vc[nv][k][j][i] - v1[nv]));
    #endif

    scrh = dt*g_maxCoolingRate/err;   
    Dts->dt_cool = MIN(Dts->dt_cool, scrh);

  /* ------------------------------------------------------
     1H. Update main solution array
     ------------------------------------------------------ */

    d->Vc[PRS][k][j][i] = prs;
    NIONS_LOOP(nv) d->Vc[nv][k][j][i] = v1[nv];

  } /* -- End DOM_LOOP(k,j,i) -- */
}
 
/* ********************************************************************* */
void Numerical_Jacobian (double *v, double **J)
/*!
 *  Compute Jacobian matrix numerically and 
 *  get eigenvector and eigenvalue decomposition
 *
 *  The purpose of this function is to detect
 *  stiffness. 
 *
 *  Note: This function is EXTREMELY time consuming.
 *
 *********************************************************************** */
{
  int  n, nv, k, l;
  double eps;
  double vpp[NVAR_COOLING], vp[NVAR_COOLING], vm[NVAR_COOLING], vmm[NVAR_COOLING];
  double Rpp[NVAR_COOLING], Rp[NVAR_COOLING], Rm[NVAR_COOLING], Rmm[NVAR_COOLING];

  eps = 1.e-5;
  n = NIONS + 1;

/* -- partial derivs with respect to pressure -- */

  for (nv = 0; nv < NVAR; nv++){
    vp[nv] = vm[nv] = v[nv];
  }
  vp[RHOE] *= 1.0 + eps;
  vm[RHOE] *= 1.0 - eps;
  
  Radiat (vp, Rp);
  Radiat (vm, Rm);

  for (k = 0; k < n - 1; k++){
    J[k][n - 1] = (Rp[k + NFLX] - Rm[k + NFLX])/(2.0*v[RHOE]*eps);
  }
  J[n - 1][n - 1] = (Rp[RHOE] - Rm[RHOE])/(2.0*v[RHOE]*eps);

/* -- partial deriv with respect to ions -- */

  for (l = 0; l < n - 1; l++){

    for (nv = 0; nv < NVAR; nv++){
      vp[nv] = vm[nv] = v[nv];
    }
    vp [l + NFLX] = v[l + NFLX] + eps;
    vm [l + NFLX] = v[l + NFLX] - eps;
                      
    vp [l + NFLX] = MIN(vp [l + NFLX], 1.0);
    vm [l + NFLX] = MAX(vm [l + NFLX], 0.0);
                         
    Radiat (vp , Rp);
    Radiat (vm , Rm);

    for (k = 0; k < n - 1; k++){
      J[k][l] = (Rp[k + NFLX] - Rm[k + NFLX])/(vp[l + NFLX] - vm[l + NFLX]); 
    }
    J[n - 1][l] = (Rp[RHOE] - Rm[RHOE])/(vp[l + NFLX] - vm[l + NFLX]);
  }
}
#else
/* --------------------------------------------------------
    			Townsend cooling
   -------------------------------------------------------- */
  
double *g_invT_tab, *g_invY_tab, *g_Y_tab, *g_L_tab, *g_T_tab;
int g_ntab;
double interp1D(double* , double* , double , char* );

/* ********************************************************************* */
void CoolingSource (const Data *d, double dt, timeStep *Dts, Grid *grid)
/*!
 * Integrate cooling and reaction source terms.
 *
 * \param [in,out]  d   pointer to Data structure
 * \param [in]     dt   the time step to be taken
 * \param [out]    Dts  pointer to the timeStep structure
 * \param [in]     grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int  nv, k, j, i, stiff, status;
  double err, scrh, min_tol = 2.e-5;
  double T0, T1, mu0, mu1, prs, tcool; 
  double Tref = 1.e9, tcool_ref= 0.;
  double v0[NVAR], v1[NVAR], k1[NVAR];
  double maxrate, dt_sub;
  int n_sub_max = 50;
  intList var_list;
  double dummy0[4], dummy1[4];
  double mue, mui;
  
  double *x  = grid->x[IDIR];
  double *y  = grid->x[JDIR];
  double *z  = grid->x[KDIR];
  double *dx = grid->dx[IDIR];
  double *dy = grid->dx[JDIR];
  double *dz = grid->dx[KDIR];
  double gCooling_x, gCooling_y, gCooling_z;
  
/* --------------------------------------------------------
    Set number and indices of the time-dependent variables
   -------------------------------------------------------- */
  
  var_list.nvar    = NIONS+1;
  var_list.indx[0] = PRS;
  for (nv = 0; nv < NIONS; nv++) var_list.indx[nv+1] = NFLX+nv;
  
/* -----------------------------------------------------
    Zero k1
   ----------------------------------------------------- */
  
  NVAR_LOOP(nv) k1[nv] = 0.0;  

  DOM_LOOP(k,j,i){  /* -- span the computational domain to find minimum tcool and global code step gets modified accordingly-- */
   NVAR_LOOP(nv) v0[nv] = v1[nv] = d->Vc[nv][k][j][i];
   
   
   mu0 = MeanMolecularWeight(v0, dummy0);
   T0  = v0[PRS]/v0[RHO]*KELVIN*mu0;
  #if EOS == IDEAL
   v0[RHOE] = v1[RHOE] = v0[PRS]/(g_gamma-1.0);
  #else
   v1[RHOE] = v0[RHOE] = InternalEnergy(v0, T0);
  #endif
   Radiat(v0, k1); //only run to populate tabulated arrays from storage if first time and get the cooling rate at every run 
   tcool = -v0[RHOE]/k1[RHOE]; //v0[RHOE]=(prs/(g_gamma-1.0)) for ideal gas and k1[RHOE]=-neniLAMBDA (See radiat.c)
   Dts->dt_cool = MIN(Dts->dt_cool, tcool);
   maxrate = GetMaxRate (v0, k1, T0);
  }
/* -----------------------------------------------------
 *     Zero reinitialization
 * ----------------------------------------------------- */
  NVAR_LOOP(nv) { k1[nv] = 0.0; v0[nv] = 0.0; v1[nv] = 0.0; }
 #ifdef PARALLEL
  double dt_cool_thisproc = Dts->dt_cool, dt_cool_allproc = 0.;
  MPI_Allreduce(&dt_cool_thisproc, &dt_cool_allproc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  Dts->dt_cool = dt_cool_allproc;
 #endif
  
/*  ----------------------------------------------------------- 
                   Begin Integration 
    -----------------------------------------------------------  */
 
/* -----------------------------------------------------
 *     Subcycling starts here
 * ----------------------------------------------------- */
  dt_sub = 1.0/(1.0/(Dts->dt_cool/n_sub_max) + 2./dt); //Subcycle limit
  int sub_steps = (int)ceil(dt/dt_sub);
  dt_sub = dt/sub_steps;
  int n_sub = 0; 
  for (n_sub = 0; n_sub<sub_steps; n_sub++ ){
    DOM_LOOP(k,j,i){  /* -- span the computational domain -- */
      /* --------------------------------------------------
         Skip integration if cell has been tagged with 
         FLAG_INTERNAL_BOUNDARY or FLAG_SPLIT_CELL 
         (only for AMR)
         -------------------------------------------------- */ 
    #if INTERNAL_BOUNDARY == YES
     if (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY) continue;
    #endif
     if (d->flag[k][j][i] & FLAG_SPLIT_CELL) continue;
     
     /* ------------------------------------------------------
     1A. Define global coordinates.
         (useful for spatial-dependent cooling
     ------------------------------------------------------ */ 

    gCooling_x = x[i];
    gCooling_y = y[j];
    gCooling_z = z[k];
     
     /* ----------------------------------------------
       Compute temperature and internal energy from
       density, pressure and concentrations.
       ---------------------------------------------- */
    
     NVAR_LOOP(nv) v0[nv] = v1[nv] = d->Vc[nv][k][j][i]; //Initialize cell fields
     prs = v0[PRS];
     mu0 = MeanMolecularWeight(v0, dummy0);
     T0  = v0[PRS]/v0[RHO]*KELVIN*mu0;
    #if EOS == IDEAL
     v0[RHOE] = v1[RHOE] = prs/(g_gamma-1.0);
    #else
     v1[RHOE] = v0[RHOE] = InternalEnergy(v0, T0);
    #endif
     
     if (T0 <= 0.0){
       print ("! CoolingSource() line 409: negative initial temperature\n");
       print (" %12.6e  %12.6e \n",v0[RHOE], v0[RHO]);
       print (" at: %f %f %f\n", gCooling_x, gCooling_y, gCooling_z);
       QUIT_PLUTO(1);
     }
     if (T0 < g_minCoolingTemp)  {
       mu1 = MeanMolecularWeight(v1, dummy1);
       prs = g_minCoolingTemp*v1[RHO]/(KELVIN*mu1);
       scrh = dt*g_maxCoolingRate/err; //****Check this out later (max change in P/ Current change in P)
       Dts->dt_cool = MIN(Dts->dt_cool, scrh); 
       d->Vc[PRS][k][j][i] = prs;
       NIONS_LOOP(nv) d->Vc[nv][k][j][i] = v1[nv]; 
       continue;
     }
     
     mu0 = MeanMolecularWeight(v0, dummy0); 
     mue = dummy0[0]; mui = dummy0[1];  
     Radiat(v0, k1); //only run to populate tabulated arrays from storage if first time and get the cooling rate at every run
     tcool = -v0[RHOE]/k1[RHOE]; //v0[RHOE]=(prs/(g_gamma-1.0)) for ideal gas and k1[RHOE]=-nenHLAMBDA (See radiat.c)
     tcool_ref = (1./(g_gamma-1))*(CONST_kB*Tref/(v0[RHO]*UNIT_DENSITY*lambda_interp(Tref)))*(mue*mui/mu0)*CONST_mp;
     tcool_ref /= (UNIT_LENGTH/UNIT_VELOCITY);
     T1 = invY_interp( Y_interp(T0) + (T0/Tref)* (lambda_interp(Tref)/lambda_interp(T0))* (dt_sub/tcool));
     mu1 = MeanMolecularWeight(v1, dummy1);
    #if EOS == IDEAL
     prs  = T1*v0[RHO]/(KELVIN*mu1);
     v1[RHOE] = prs/(g_gamma - 1.0);
    #else
     v1[RHOE] = InternalEnergy(v1, T1);
    #endif
     /* -- pressure must be positive -- */
     if (prs < 0.0) prs = g_smallPressure;
     /* ------------------------------------------------------
      *       As a safety measure, constrain the ions to
      *               lie in [0,1] and pressure to be positive.
      * ------------------------------------------------------ */
     NIONS_LOOP(nv){
       v1[nv] = MAX(v1[nv], 0.0);
       v1[nv] = MIN(v1[nv], 1.0);
     }
     /* ------------------------------------------------------
      *       If gas is initially or cooled below lower threshold, redefine
      *               pressure so that T = g_minCoolingTemp.
      * ------------------------------------------------------ */
     if (T1 < g_minCoolingTemp || T0 < g_minCoolingTemp)  prs = g_minCoolingTemp*v1[RHO]/(KELVIN*mu1);
     
     err = fabs(prs/d->Vc[PRS][k][j][i] - 1.0);
    #if COOLING == MINEq
     for (nv = NFLX; nv < NFLX + NIONS - Fe_IONS; nv++) 
    #else
     NIONS_LOOP(nv)  err = MAX(err, fabs(d->Vc[nv][k][j][i] - v1[nv]));
    #endif
     err = MAX(err, fabs(d->Vc[nv][k][j][i] - v1[nv]));
     
     scrh = dt*g_maxCoolingRate/err; //****Check this out later (max change in P/ Current change in P)
     
     Dts->dt_cool = MIN(Dts->dt_cool, scrh);
     
     /* ------------------------------------------------------
      *      1H. Update main solution array
      * ------------------------------------------------------ */
     d->Vc[PRS][k][j][i] = prs;
     NIONS_LOOP(nv) d->Vc[nv][k][j][i] = v1[nv];
    } /* -- End DOM_LOOP(k,j,i) -- */
  //print("Temp %e     step %d     time %e     tcool %e\n",T1,n_sub,g_time,Dts->dt_cool);
  }   /* -- End Subcycle Loop   -- */
}

double lambda_interp(double temperature) {
  return interp1D(g_T_tab, g_L_tab, temperature, "lambda_interp");
}

double Y_interp(double temperature) {
  return interp1D(g_T_tab, g_Y_tab, temperature, "Y_interp");
}

double invY_interp(double townY) {
  return interp1D(g_invY_tab, g_invT_tab, townY, "invY_interp");
}

double interp1D(double* x_data, double* y_data, double x_request, char* msg) {

 /* ----------------------------------------------
  *         Table lookup by binary search used for interpolating y as a function of x. 
  *        x is assumed to be arranged in ascending order
  *  ---------------------------------------------- */
  int ntab = g_ntab; // number of entries maybe (int)(sizeof(x_data)/sizeof(x_data[0])) will work
  int klo = 0;
  int khi = ntab - 1;
  int kmid;
  double xmid, dx, scrh;

  if (x_request > x_data[khi] || x_request < x_data[klo]){
    print ("Called from %s\n",msg);
    print (" ! Requested value out of range: %12.6e\n",x_request);
    QUIT_PLUTO(1);
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    xmid = x_data[kmid];
    if (x_request <= xmid){
      khi = kmid;
    }else if (x_request > xmid){
      klo = kmid;
    }
  }

/* -----------------------------------------------
 *     Compute r.h.s
 * ----------------------------------------------- */

  dx       = x_data[khi] - x_data[klo];
  scrh     = y_data[klo]*(x_data[khi] - x_request)/dx + y_data[khi]*(x_request - x_data[klo])/dx;

  return scrh;
}
#endif
