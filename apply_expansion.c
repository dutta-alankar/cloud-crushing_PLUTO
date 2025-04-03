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
  
  /* Expand the grid */
  KTOT_LOOP(k){
    grid->xr_glob[KDIR][k] += (0.5*scale*grid->dx_glob[KDIR][k]);
    grid->xl_glob[KDIR][k] -= (0.5*scale*grid->dx_glob[KDIR][k]);
    grid->x_glob[KDIR][k] = 0.5*(grid->xl_glob[KDIR][k]+grid->xr_glob[KDIR][k]);
    grid->dx_glob[KDIR][k] *= scale;
  }
  JTOT_LOOP(j){
    grid->xr_glob[JDIR][j] += (0.5*scale*grid->dx_glob[KDIR][j]);
    grid->xl_glob[JDIR][j] -= (0.5*scale*grid->dx_glob[KDIR][j]);
    grid->x_glob[JDIR][j] = 0.5*(grid->xl_glob[KDIR][j]+grid->xr_glob[KDIR][j]);
    grid->dx_glob[JDIR][j] *= scale;
  }
  grid->xbeg_glob[KDIR] = g_domBeg[KDIR] = grid->xl_glob[KDIR][grid->gbeg[KDIR]];
  grid->xbeg_glob[JDIR] = g_domBeg[JDIR] = grid->xl_glob[JDIR][grid->gbeg[JDIR]];
  grid->xend_glob[KDIR] = g_domEnd[KDIR] = grid->xr_glob[KDIR][grid->gend[KDIR]];
  grid->xend_glob[JDIR] = g_domEnd[JDIR] = grid->xr_glob[JDIR][grid->gend[JDIR]];
  SetGeometry(grid);
  
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
