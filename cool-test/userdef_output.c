#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k; 
  double *x1, *x2, *x3;
  double ***T;

  /* get grid pointers */
  x1  = grid->x[IDIR];   x2 = grid->x[JDIR];   x3 = grid->x[KDIR];
 
  #if COOLING==NO || COOLING==TABULATED || COOLING==TOWNSEND
  double dummy[4];
  double mu = MeanMolecularWeight((double*)d->Vc, dummy);
  #else
  double mu = MeanMolecularWeight((double*)d->Vc);
  #endif
 
  T    = GetUserVar("Temp");

  TOT_LOOP(k,j,i) {
     T[k][j][i] = (d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i]) * pow(UNIT_VELOCITY,2) * (mu * CONST_mp)/CONST_kB;
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#if PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}
