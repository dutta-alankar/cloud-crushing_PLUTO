#include "pluto.h"
#include "local_pluto.h"

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
  double *dx, *dy, *dz;

  double rho, dV;

  double ***celldV;
  double *dx1, *dx2, *dx3;
  double ***T, ***n, ***p, ***mach;

  /* get grid pointers */
  x1  = grid->x[IDIR];   x2 = grid->x[JDIR];   x3 = grid->x[KDIR];
  dx1 = grid->dx[IDIR]; dx2 = grid->dx[JDIR]; dx3 = grid->dx[KDIR];
 
  #if COOLING==NO || COOLING==TABULATED || COOLING==TOWNSEND
  double dummy[4];
  double mu = MeanMolecularWeight((double*)d->Vc, dummy);
  #else
  double mu = MeanMolecularWeight((double*)d->Vc);
  #endif
 
  T    = GetUserVar("Temp");
  n    = GetUserVar("ndens");
  p    = GetUserVar("PbykB");
  mach = GetUserVar("mach");
  celldV = GetUserVar("cellvol");

  TOT_LOOP(k,j,i) {
     rho   = d->Vc[RHO][k][j][i];
     dV    = grid->dV[k][j][i];
     T[k][j][i] = (d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i]) * pow(UNIT_VELOCITY,2) * (mu * CONST_mp)/CONST_kB;
     n[k][j][i] = d->Vc[RHO][k][j][i] * UNIT_DENSITY / (mu * CONST_mp);
     p[k][j][i] = d->Vc[PRS][k][j][i] * UNIT_DENSITY * pow(UNIT_VELOCITY,2.) / CONST_kB;
     mach[k][j][i] = sqrt( DIM_EXPAND(d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i],
                                  + d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i],
                                  + d->Vc[VX3][k][j][i]*d->Vc[VX3][k][j][i]) )/sqrt(g_gamma*(d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]));
     celldV[k][j][i] = grid->dV[k][j][i];                            
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
