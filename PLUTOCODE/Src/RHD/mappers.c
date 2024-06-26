/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the RHD equations.
  
  The ConsToPrim() converts an array of conservative quantities to 
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total 
  energy unless zone has been tagged with FLAG_ENTROPY.
  In this case we recover pressure from conserved entropy:
  
      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)
  
  \author A. Mignone (mignone@to.infn.it)
  \date   May 3, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define  ENERGY_SOLVE   1
#define  ENTROPY_SOLVE  2
#define  PRESSURE_FIX   3
/* ********************************************************************* */
void PrimToCons  (double *uprim[], double *ucons[],
                 int ibeg, int iend)
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
  int     nv, ii;
  double  rhoh_g2, scrh, g;
  double  beta_fix = 0.9999;
  double  *u, *v;
  static double  *h;

  if (h == NULL){
    h = ARRAY_1D(NMAX_POINT, double);
  }

  Enthalpy (uprim, h, ibeg, iend);

  for (ii = ibeg; ii <= iend; ii++) {
   
    u = ucons[ii];
    v = uprim[ii];

    g = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];

    if (g >= 1.0){
      WARNING( 
        printLog ("! u^2 > 1 (%f) in PrimToCons\n", scrh);
        Where (ii, NULL);
      )

      g = beta_fix/sqrt(g);
      v[VX1] *= g;
      v[VX2] *= g;
      v[VX3] *= g;

      g = beta_fix*beta_fix;
    }
    g    = 1.0/(1.0 - g);
    scrh = rhoh_g2 = v[RHO]*h[ii]*g;
    g    = sqrt(g);

    u[RHO] = v[RHO]*g;
    u[MX1] = scrh*v[VX1];
    u[MX2] = scrh*v[VX2];
    u[MX3] = scrh*v[VX3];

    ucons[ii][ENG] = rhoh_g2 - v[PRS];

    #if NSCL > 0    
    NSCL_LOOP(nv) u[nv] = v[nv]*u[RHO];
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
int ConsToPrim (double **ucons, double **uprim, int ibeg, int iend, 
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
  int     nv, i, err, ifail;
  char    sol;
  double  *u, *v;

  ifail = 0;
  for (i = ibeg; i <= iend; i++) {

    u = ucons[i];
    v = uprim[i];

    sol = ENERGY_SOLVE;
    #if ENTROPY_SWITCH
    if (flag[i] & FLAG_ENTROPY) sol = ENTROPY_SOLVE;
    #endif

  /* --------------------------------------------
     1. Check density and rad energy positivity 
     -------------------------------------------- */
  
    if (u[RHO] < 0.0) {
      printLog("! ConsToPrim(): negative density (%8.2e), ", u[RHO]);
      Where (i, NULL);
      u[RHO]   = g_smallDensity;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }
    
    #if RADIATION
    if (u[ENR] < 0.0) {
      WARNING(
        print("! ConsToPrim(): negative radiation energy (%8.2e), ", u[ENR]);
        Where (i, NULL);
      )
      u[ENR]   = RADIATION_MIN_ERAD;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }
    #endif

  /* --------------------------------------------
     4a. Obtain pressure by inverting 
         the entropy equation.
     -------------------------------------------- */

    #if ENTROPY_SWITCH
    if (sol == ENTROPY_SOLVE) {
      err = RHD_EntropySolve(u, v);
      if (err) {
        WARNING(
          printLog ("! ConsToPrim() <- RHD_EntropySolve() failed, ");
          printLog ("err code = %d. Trying pressure fix... \n", err);
          Where(i,NULL);
        )
        sol = PRESSURE_FIX;
      }
    }  /* end if entropy_solve */
    #endif
 
  /* --------------------------------------------
     4b. Obtain pressure by inverting 
         the energy equation.
     -------------------------------------------- */

    if (sol == ENERGY_SOLVE){
      err = RHD_EnergySolve(u, v);
      if (err){
        WARNING(
          printLog ("! ConsToPrim(): RHD_EnergySolve() failed,");
          printLog (" err code = %d; ", err);
          Where(i,NULL);
        )
        sol = PRESSURE_FIX;
      }
    } /* end if energy_solve */

  /* --------------------------------------------
     4c. Apply pressure fix if necessary
     -------------------------------------------- */

    if (sol == PRESSURE_FIX){
      err = RHD_PressureFix(u, v);
      if (err){
        printLog ("! ConsToPrim() <- RHD_PressureFix() failed, ");
        printLog ("err code = %d; \n", err);
        Where(i,NULL);
        ShowState (u,0);          
        QUIT_PLUTO(1);  /* No more options  (!!) */
      }
      
    /* -- Pressure has changed, redefine total energy & entropy -- */
      
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    } /* end if pressure_fix */

  /* --------------------------------------------
     5. Complete conversion 
     -------------------------------------------- */

    #if NSCL > 0     
    NSCL_LOOP(nv)  v[nv] = u[nv]/u[RHO];
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
