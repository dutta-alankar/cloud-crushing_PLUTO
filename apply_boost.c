#include "pluto.h"
#include "local_pluto.h"

void ApplyBoost (const Data *d, double dt, timeStep *Dts, Grid *grid)
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

  x1  = grid->x[IDIR];   x2 = grid->x[JDIR];   x3 = grid->x[KDIR];
  dx1 = grid->dx[IDIR]; dx2 = grid->dx[JDIR]; dx3 = grid->dx[KDIR];
  
  /* --- variables used --- */
  double rho, vx1;
  double vx1_cl = 0.;
  double mass = 0.;
  double tracerLeftEdge = 2.0*grid->xend_glob[IDIR];
  
  /* ----- Main Loop ----- */
  DOM_LOOP(k,j,i){
    dV  = grid->dV[k][j][i];
    rho = d->Vc[RHO][k][j][i];
    vx1 = d->Vc[VX1][k][j][i];  // z-velocity 

    if ((grid->xl[IDIR][i] < tracerLeftEdge) && (d->Vc[TRC][k][j][i] >= 1e-4))
        tracerLeftEdge = grid->xl[IDIR][i];

    vx1_cl += vx1 * rho * d->Vc[TRC][k][j][i] * dV;
    mass   += rho * d->Vc[TRC][k][j][i] * dV;
  }
  
  #ifdef PARALLEL
  int count = 2;
  int tmp = 0;
  double sendArray[count], recvArray[count];

  sendArray[tmp++] = vx1_cl;  sendArray[tmp++] = mass; 
    
  // MPI_Reduce(sendArray, recvArray, count, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Allreduce(sendArray, recvArray, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  tmp = 0;
  vx1_cl = recvArray[tmp++];   mass   = recvArray[tmp++];
  #endif // end of PARALLEL

  vx1_cl /= mass; /* com velocity */
  
  #ifdef PARALLEL
  double dummy;
  MPI_Allreduce(&tracerLeftEdge, &dummy, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  tracerLeftEdge = dummy;
  #endif
  g_tracerLeftEdge =  tracerLeftEdge;
  
  if ( (fabs(tracerLeftEdge-grid->xbeg_glob[IDIR])<8.0) && (vx1_cl>0) ) {
    g_vcloud += vx1_cl; // frame velocity wrt lab

    TOT_LOOP(k,j,i){
      d->Vc[VX1][k][j][i] -= vx1_cl;
    }
  
    #ifdef PARALLEL
    /* This call communicates the field update across processor */
    Boundary (d, ALL_DIR, grid);
    MPI_Barrier (MPI_COMM_WORLD);
    #endif
    RBox dom_box;
    RBoxDefine (0, NX1_TOT-1, 0,  NX2_TOT-1, 0, NX3_TOT-1, CENTER, &dom_box);
    PrimToCons3D (d->Vc, d->Uc, &dom_box);
  }
  return ;
}
