#include "pluto.h"
#include "local_pluto.h"
#include <math.h>
extern double *g_invT_tab, *g_invY_tab, *g_Y_tab, *g_L_tab, *g_T_tab;
extern int g_ntab;

/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
{
  //int    klo, khi, kmid;
  static int ntab;
  double  mu, muH, mue, mui, T, lam, ne, ni, nH, n, prs; 
  double E_cost;
  double dummy[4];
  
  FILE *fcool;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

  if (g_T_tab == NULL){ //This line ensures reading is done only once per run
    print (" > Reading table %s from disk...\n","cooltable_townsend.dat");
    fcool = fopen("./cooltable_townsend.dat","r");
    if (fcool == NULL){
      print ("! Radiat: %s could not be found.\n","cooltable_townsend.dat");
      QUIT_PLUTO(1);
    }
    g_L_tab =    ARRAY_1D(20000, double);
    g_T_tab =    ARRAY_1D(20000, double);
    g_Y_tab =    ARRAY_1D(20000, double); 
    g_invT_tab = ARRAY_1D(20000, double);
    g_invY_tab = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf %lf %lf %lf %lf\n", g_T_tab + ntab, 
                                       g_L_tab + ntab, g_Y_tab + ntab, g_invY_tab + ntab, g_invT_tab + ntab)!=EOF) { 
	  ntab++;
    }
    g_ntab = ntab; //number of rows in cooling data
  }
  E_cost = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH;
  //print("mark %e\n",g_invY_tab[0]);
/* ---------------------------------------------
            Get pressure and temperature 
   --------------------------------------------- */

  prs = v[PRS]; //v[RHOE]*(g_gamma-1.0);
  if (prs < 0.0) {
    prs     = g_smallPressure;
    v[RHOE] = prs/(g_gamma - 1.0); //internal thermal energy
  }

  mu  = MeanMolecularWeight(v, dummy);
  T = prs/v[RHO]*KELVIN*mu;
  mue = dummy[0]; mui = dummy[1]; muH = dummy[2];
  /*
  if (T != T){
    printf (" ! Nan found in radiat \n");
    printf (" ! rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
    QUIT_PLUTO(1);
  }
  
  if (T < g_minCoolingTemp) { 
    rhs[RHOE] = 0.0;
    return;
  }
  */

/* -----------------------------------------------
    Compute r.h.s
   ----------------------------------------------- */

  lam      = lambda_interp(T);
  
  ne       = v[RHO]*UNIT_DENSITY/(CONST_mp*mue); //CGS
  ni       = v[RHO]*UNIT_DENSITY/(CONST_mp*mui); //CGS
  nH       = v[RHO]*UNIT_DENSITY/(CONST_mp*muH); //CGS
  n        = v[RHO]*UNIT_DENSITY/(CONST_mp*mu); //CGS
  rhs[RHOE] = -nH*nH*lam; //ne ni LAMBDA in CGS
  rhs[RHOE] /=  E_cost; ///UPDATE Cooling loss in code units

/* ----------------------------------------------
 *     Temperature cutoff
 * ---------------------------------------------- */

  rhs[RHOE] *= 1.0 - 1.0/cosh( pow( T/g_minCoolingTemp, 12)); //ram down Lambda
}
