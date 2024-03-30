/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Return a pointer to a Riemann solver function.

  \author A. Mignone (mignone@to.infn.it)
  \date   Sep 17, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
Riemann_Solver *SetSolver (const char *solver)
/*!
 *  Depending on the choice of the Riemann solver specified in
 *  pluto.ini, return a pointer to the corresponding Riemann solver
 *  function
 *
 *********************************************************************** */
{
  
/* ------------------------------------------------------
       Set Pointers for SOLVERS 
   ------------------------------------------------------ */

  if (!strcmp(solver, "tvdlf"))         return (&LF_Solver);
  else if (!strcmp(solver, "hlle") || 
           !strcmp(solver, "hll"))      return (&HLL_Solver);
  else if (!strcmp(solver, "hllc"))     return (&HLLC_Solver);
  else if (!strcmp(solver, "hlld"))     return (&HLLD_Solver);
  else if (!strcmp(solver, "gforce"))   return (&GFORCE_Solver);
/*
  else if (!strcmp(solver, "gmusta1"))  return (&GMUSTA1_Solver);
  else if (!strcmp(solver, "musta"))      return (&MUSTA_Solver);
  else if (!strcmp(solver, "rusanov_dw")) return (&RusanovDW_Solver);
*/

  printLog ("\n ! SetSolver: '%s' is not available.\n", solver);
  QUIT_PLUTO(1);

}
