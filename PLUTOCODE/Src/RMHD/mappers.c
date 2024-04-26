/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the RMHD equations.
  
  The ConsToPrim() converts an array of conservative quantities to 
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total 
  energy using the algorithm outlined in
  - "Equation of sweep in relativistic magnetohydrodynamics: variable versus
     constant adiabatic index"\n
     Mignone \& Mc Kinney, MNRAS (2007) 378, 1118.

  However, if the zone has been tagged with FLAG_ENTROPY, primitive
  variables are recovered by using the conserved entropy rather than
  total energy.
  
  In other words:
  
      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)
  
  If the inversion scheme fails and p cannot be obtained the
  RMHD_PressureFix() function is used.

  \author A. Mignone (mignone@to.infn.it)
  \date   July 1, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void PrimToCons (double **uprim, double **ucons, int beg, int end)
/*!
 * Convert primitive variables to conservative variables. 
 *
 * \param [in]  uprim array of primitive variables
 * \param [out] ucons array of conservative variables
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final index of computation
 *
 *********************************************************************** */
{
  int   i, nv;
  double  vel2, usq, vB, Bmag2;
  double  vx1, vx2, vx3;
  double  g, g2, wt;
  double  *u, *v;
  static double *h;
  #if EOS == IDEAL
  double gmmr = g_gamma/(g_gamma - 1.0);
  #endif

  if (h == NULL) h = ARRAY_1D(NMAX_POINT, double);

  Enthalpy(uprim, h, beg, end);

  for (i = beg; i <= end; i++) {

    v = uprim[i];
    u = ucons[i];

    vx1  = v[VX1];
    vx2  = v[VX2];
    vx3  = v[VX3];
    vel2 = vx1*vx1 + vx2*vx2 + vx3*vx3;
    g2 = 1.0/(1.0 - vel2);
    g  = sqrt(g2);

    vB    = vx1*v[BX1]    + vx2*v[BX2]    + vx3*v[BX3];
    Bmag2 = v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3];
    wt    = v[RHO]*h[i]*g2 + Bmag2;

  /* -------------------------------------------------------
       Convert from primitive (v) to conservative (u)   
     ------------------------------------------------------- */

    u[RHO] = g*v[RHO];
    u[MX1] = wt*vx1 - vB*v[BX1];
    u[MX2] = wt*vx2 - vB*v[BX2];
    u[MX3] = wt*vx3 - vB*v[BX3];

    u[BX1] = v[BX1];
    u[BX2] = v[BX2];
    u[BX3] = v[BX3];

    #if RMHD_REDUCED_ENERGY == YES
    #if EOS == IDEAL
    u[ENG]  =   v[PRS]*(g2*gmmr - 1.0) + u[RHO]*g2*vel2/(g + 1.0)
              + 0.5*(Bmag2*(1.0 + vel2) - vB*vB);
    #elif EOS == TAUB
    wt     =   v[PRS]/v[RHO];
    u[ENG] =   v[PRS]*(g2*2.5 - 1.0) 
             + u[RHO]*g2*(2.25*wt*wt + vel2)/(g*sqrt(1.0 + 2.25*wt*wt) + 1.0)
             + 0.5*(Bmag2*(1.0 + vel2) - vB*vB);
    #endif
    #else
    u[ENG]  = v[RHO]*h[i]*g2 - v[PRS] + 0.5*(Bmag2*(1.0 + vel2) - vB*vB);
    #endif  /* RMHD_REDUCED_ENERGY */

    #if NSCL > 0
    NSCL_LOOP(nv) u[nv] = u[RHO]*v[nv];
    #endif    
    
    #ifdef GLM_MHD
     u[PSI_GLM] = v[PSI_GLM]; 
    #endif

	  #if RADIATION
		u[ENR] = v[ENR];
		u[FR1] = v[FR1]; 
    u[FR2] = v[FR2]; 
    u[FR3] = v[FR3];
	  #endif
  }
}

/* ********************************************************************* */
int ConsToPrim (double **ucons, double **uprim, int beg, int end, 
                unsigned char *flag)
/*!
 * Convert from conservative to primitive variables.
 *
 * \param [in]  ucons      array of conservative variables
 * \param [out] uprim      array of primitive variables
 * \param [in]  beg        starting index of computation
 * \param [in]  end        final index of computation
 * \param [in,out] flag    array of flags tagging, in input, zones
 *                         where entropy must be used to recover pressure
 *                         and, on output, zones where conversion was
 *                         not successful.
 * 
 * \return Return 0 if conversion was successful in all zones in the 
 *         range [ibeg,iend].
 *         Return 1 if one or more zones could not be converted correctly
 *         and either pressure, density or energy took non-physical values. 
 *
 *********************************************************************** */
{
  int    i, nv, err, ifail;
  int    use_entropy, use_energy=1;
  double *u, *v, scrh, w_1;
  Map_param par;

  ifail = 0;
  for (i = beg; i <= end; i++) {

    u = ucons[i];
    v = uprim[i];

  /* --------------------------------------------
     1, Define the input parameters of the
        parameter structure
     -------------------------------------------- */

    par.D  = u[RHO];
    par.E  = u[ENG];
    par.S  = u[MX1]*u[BX1] + u[MX2]*u[BX2] + u[MX3]*u[BX3];
    par.m2 = u[MX1]*u[MX1] + u[MX2]*u[MX2] + u[MX3]*u[MX3];
    par.B2 = u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3]; 
    par.S2 = par.S*par.S;

  /* --------------------------------------------
     2. Check density and energy positivity 
     -------------------------------------------- */
  
    if (u[RHO] < 0.0) {
      printLog("! ConsToPrim(): negative density (%8.2e), ", u[RHO]);
      Where (i, NULL);
      u[RHO]   = g_smallDensity;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }

    if (u[ENG] < 0.0) {
      WARNING(
        printLog("! ConsToPrim(): negative energy (%8.2e), ", u[ENG]);
        Where (i, NULL);
      )
      u[ENG]   = 1.e-5;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }

    #if RADIATION
    if (u[ENR] < 0.0) {
      WARNING(
        print("! ConsToPrim: negative radiation energy (%8.2e), ", u[ENR]);
        Where (i, NULL);
      )			
      u[ENR]   = RADIATION_MIN_ERAD;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }
    #endif

  /* --------------------------------------------
     3. Attempt to recover pressure and velocity
        from energy or entropy.
        If an error occurs, use the
        RMHD_PressureFix() function
     -------------------------------------------- */

    #if ENTROPY_SWITCH
    use_entropy = (flag[i] & FLAG_ENTROPY);
    use_energy  = !use_entropy;
    par.sigma_c = u[ENTR];
    if (use_entropy) {
      err = RMHD_EntropySolve(&par);      
      if (err) {
        WARNING(Where (i, NULL);)
        err = RMHD_PressureFix(&par);
        if (err){
          Where(i,NULL);
        }
        u[ENTR]  = par.sigma_c;
        flag[i] |= FLAG_CONS2PRIM_FAIL;
        ifail    = 1;
      }
      u[ENG] = par.E;  /* Redefine energy */
    } 
    #endif
 
    if (use_energy){
      err = RMHD_EnergySolve(&par);

      if (err){
        WARNING(
          printLog ("! ConsToPrim(): RMHD_EnergySolve() failed,");
          printLog (" err code = %d;", err);
          Where(i,NULL);
        )

        err = RMHD_PressureFix(&par);
        if (err){
          printLog ("! ConsToPrim(): RMHD_PressureFix() failed,");
          printLog (" err code = %d;", err);
          Where(i,NULL);
//          QUIT_PLUTO(1);
        }
        u[ENG]   = par.E;
        flag[i] |= FLAG_CONS2PRIM_FAIL;
        ifail    = 1;
      }
      #if ENTROPY_SWITCH     
      u[ENTR] = par.sigma_c;  /* Redefine entropy */
      #endif      
    }

  /* --------------------------------------------
     4. W, p and lor have been found.
        Now complete conversion
     -------------------------------------------- */

    v[RHO] = u[RHO]/par.lor;
    v[PRS] = par.prs;

    w_1  = 1.0/(par.W + par.B2);
    scrh = par.S/par.W;
    v[VX1] = w_1*(u[MX1] + scrh*u[BX1]);
    v[VX2] = w_1*(u[MX2] + scrh*u[BX2]);
    v[VX3] = w_1*(u[MX3] + scrh*u[BX3]);

    scrh = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
    if (scrh >= 1.0){
      printLog ("! ConsToPrim(): v^2 = %f > 1  (p = %12.6e); ", scrh, par.prs);
      Where (i, NULL);
      printLog ("!               Flag_Entropy = %d\n", (flag[i] & FLAG_ENTROPY));
      flag[i] |= FLAG_CONS2PRIM_FAIL;

//      QUIT_PLUTO(1);
    }

    v[BX1] = u[BX1];
    v[BX2] = u[BX2];
    v[BX3] = u[BX3];

    #if NSCL > 0 
    NSCL_LOOP(nv) v[nv] = u[nv]/u[RHO];
    #endif

    #ifdef GLM_MHD
    v[PSI_GLM] = u[PSI_GLM]; 
    #endif

	  #if RADIATION
		v[ENR] = u[ENR];
		v[FR1] = u[FR1]; 
    v[FR2] = u[FR2]; 
    v[FR3] = u[FR3];
	  #endif

  }
  return ifail;
}
