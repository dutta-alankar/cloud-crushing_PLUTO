#include "pluto.h"
#include "local_pluto.h"

void ApplyExpansion (const Data *d, double dt, timeStep *Dts, Grid *grid)
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
  int nv, k, j, i;
  double *x1, *x2, *x3;
  double *dx1, *dx2, *dx3;
  double dV;
  static double scale_new = 1.0; 
  static double scale_old = 1.0;
  double scale;

  x1  = grid->x[IDIR];   x2 = grid->x[JDIR];   x3 = grid->x[KDIR];
  dx1 = grid->dx[IDIR]; dx2 = grid->dx[JDIR]; dx3 = grid->dx[KDIR];

  double x_offset = g_inputParam[XOFFSET];
  scale_new = g_dist_lab/x_offset;
  scale = scale_new/scale_old;
  scale_old = scale_new;

  double rho_pow = -2.0;
  double prs_pow = -2.0*g_gamma;
  double vpr_pow = -1.0;
  /* ----- Main Loop ----- */
  TOT_LOOP(k,j,i){
    d->Vc[RHO][k][j][i] *= pow(scale,  rho_pow);
    d->Vc[PRS][k][j][i] *= pow(scale,  prs_pow);
    d->Vc[VX2][k][j][i] *= pow(scale,  vpr_pow);
    d->Vc[VX3][k][j][i] *= pow(scale,  vpr_pow);
    d->Vc[TRC][k][j][i] *= pow(scale, -rho_pow);
  }
  
  #ifdef PARALLEL
  /* This call communicates the field update across processor */
  Boundary (d, ALL_DIR, grid);
  MPI_Barrier (MPI_COMM_WORLD);
  #endif
  RBox dom_box;
  RBoxDefine (0, NX1_TOT-1, 0,  NX2_TOT-1, 0, NX3_TOT-1, CENTER, &dom_box);
  PrimToCons3D (d->Vc, d->Uc, &dom_box);
  return ;
}
