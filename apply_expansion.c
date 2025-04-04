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
    // d->Vc[TRC][k][j][i] *= pow(scale, -rho_pow);
  }
  
  /* Expand the grid */
  for(k=0; k<grid->np_tot_glob[KDIR]; k++){
    grid->xr_glob[KDIR][k] *= scale;
    grid->xl_glob[KDIR][k] *= scale;
    grid->x_glob[KDIR][k]   = (0.5*(grid->xl_glob[KDIR][k]+grid->xr_glob[KDIR][k]));
    grid->dx_glob[KDIR][k] *= scale;
  }
  for(j=0; j<grid->np_tot_glob[JDIR]; j++){
    grid->xr_glob[JDIR][j] *= scale;
    grid->xl_glob[JDIR][j] *= scale;
    grid->x_glob[JDIR][j]   = (0.5*(grid->xl_glob[JDIR][j]+grid->xr_glob[JDIR][j]));
    grid->dx_glob[JDIR][j] *= scale;
  }
  /* Re-assign the global beg and end of the domain */
  grid->xbeg_glob[KDIR] = g_domBeg[KDIR] = grid->xl_glob[KDIR][grid->gbeg[KDIR]];
  grid->xbeg_glob[JDIR] = g_domBeg[JDIR] = grid->xl_glob[JDIR][grid->gbeg[JDIR]];
  grid->xend_glob[KDIR] = g_domEnd[KDIR] = grid->xr_glob[KDIR][grid->gend[KDIR]];
  grid->xend_glob[JDIR] = g_domEnd[JDIR] = grid->xr_glob[JDIR][grid->gend[JDIR]];
  /* Re-calculate areas and volume for every cell */
  SetGeometry(grid);

  RBox tot_box;
  RBoxDefine (0, NX1_TOT-1, 0,  NX2_TOT-1, 0, NX3_TOT-1, CENTER, &tot_box);
  PrimToCons3D (d->Vc, d->Uc, &tot_box);
  
  #ifdef PARALLEL
  /* This call communicates the field update across processors */
  Boundary (d, ALL_DIR, grid);
  MPI_Barrier (MPI_COMM_WORLD);
  #endif
  return ;
}
