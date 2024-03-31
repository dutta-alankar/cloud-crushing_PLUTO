/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "local_pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  g_gamma      = 5/3.;
  g_tracerLeftEdge = -1.;
  
  /* set pointer shortcuts */
  int i, j, k;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  double oth_mu[4];
  double mu = MeanMolecularWeight((double*)d->Vc, oth_mu);
  printLog("> Mean molecular weights: \n"); 
  printLog("  mu = %.3f\tmue = %.3f\tmui = %.3f\n\n",mu, oth_mu[0], oth_mu[1]);

  double nmin     = 1e-6;
  g_smallDensity  = (nmin*((CONST_mp*mu)/UNIT_DENSITY));
  g_smallPressure = (nmin*g_minCoolingTemp*CONST_kB)/(UNIT_DENSITY*pow(UNIT_VELOCITY,2));
  
  /* wind properties in code units */
  double rhoWind = 1.0; 
  double mach    = g_inputParam[MACH];
  double pWind   = 1.0/(g_gamma*pow(mach, 2.));
  double vWind   = 1.0;  

  /* cloud properties in code units */
  double chi      = g_inputParam[CHI];
  double eta      = g_inputParam[ETA];
  double x_offset = g_inputParam[XOFFSET];

  /* ---- Main Loop ---- */
  TOT_LOOP(k,j,i) {
    double distance = sqrt(pow(x1[i]-x_offset, 2.) + pow(x2[j], 2.) + pow(x3[k], 2.));
    d->Vc[RHO][k][j][i] = ((distance <= 1.0)? chi : 1.) * rhoWind; 
    d->Vc[PRS][k][j][i] = ((distance <= 1.0)? (chi/eta) : 1.) * pWind;
    DIM_EXPAND(
      d->Vc[VX1][k][j][i] = (distance > 1.0)? vWind : 0.;,
      d->Vc[VX2][k][j][i] = 0.;,
      d->Vc[VX3][k][j][i] = 0.;
    )
    d->Vc[TRC][k][j][i] = (distance <= 1.0)? 1.0 : 0.;
  } // end of TOT_LOOP
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
  int i, j, k;
  double *x1, *x2, *x3;
  double *dx1, *dx2, *dx3;
  #ifdef PARALLEL
  int count = 128;
  double sendArray[count], recvArray[count];
  #endif
  
  /* set grid pointers */
  x1  = grid->x[IDIR];   x2 = grid->x[JDIR];   x3 = grid->x[KDIR];
  dx1 = grid->dx[IDIR]; dx2 = grid->dx[JDIR]; dx3 = grid->dx[KDIR];

  /* variables used */
  double rho, prs, vx1, vx2, vx3, Temp, dV;

  double tot_mass = 0.,  tot_mom1 = 0.,  tot_mom2 = 0.,  tot_mom3 = 0.;
  double tot_TE = 0., tot_vol=0.;
  double vx1_avg = 0.,   vx2_avg = 0.,   vx3_avg = 0.,   v_rms = 0.;

  double dummy[4];
  double mu = MeanMolecularWeight((double*)d->Vc, dummy);
  
  /* wind properties in code units */
  double rhoWind = 1.0; 
  double mach    = g_inputParam[MACH];
  double pWind   = 1.0/(g_gamma*pow(mach, 2.));
  double vWind   = 1.0;  

  /* cloud properties in code units */
  double chi      = g_inputParam[CHI];
  double x_offset = g_inputParam[XOFFSET];

  /* calculating global averages first */
  DOM_LOOP(k,j,i) {
    dV  = grid->dV[k][j][i];
    rho = d->Vc[RHO][k][j][i];
    prs = d->Vc[PRS][k][j][i];
    vx1 = d->Vc[VX1][k][j][i];  
    vx2 = d->Vc[VX2][k][j][i]; 
    vx3 = d->Vc[VX3][k][j][i];

    tot_vol  += dV;
    vx1_avg += vx1 * (1-d->Vc[TRC][k][j][i])* dV; // average wind velocity
    vx2_avg += vx2 * (1-d->Vc[TRC][k][j][i])* dV;
    vx3_avg += vx3 * (1-d->Vc[TRC][k][j][i])* dV;
    
    tot_TE   += prs * (1-d->Vc[TRC][k][j][i])* dV / (g_gamma-1.0);
  }
  /* ------ Parallel data reduction ------ */
  #ifdef PARALLEL
  int tmp = 0;
  count = 5;
  sendArray[tmp++] = tot_vol;
  sendArray[tmp++] = vx1_avg;   sendArray[tmp++] = vx2_avg;    sendArray[tmp++] = vx3_avg;
  sendArray[tmp++] = tot_TE;

  MPI_Allreduce(sendArray, recvArray, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  tmp = 0;
  tot_vol = recvArray[tmp++];
  vx1_avg = recvArray[tmp++];    vx2_avg = recvArray[tmp++];   vx3_avg = recvArray[tmp++];  
  tot_mom1 = recvArray[tmp++];   tot_mom2 = recvArray[tmp++];  tot_mom3 = recvArray[tmp++];
  tot_TE  = recvArray[tmp++];

  MPI_Barrier(MPI_COMM_WORLD);
  #endif
 
  vx1_avg /= tot_vol;    vx2_avg /= tot_vol;    vx3_avg /= tot_vol;

  double cloud_mass = 0., cloud_vol = 0., trc_mass = 0.;
  double cmx1_cl = 0., cmx2_cl = 0., cmx3_cl = 0.;
  double vx1_cl = 0., vx2_cl = 0., vx3_cl = 0.;

  tot_vol = 0.;
  DOM_LOOP(k,j,i) {
    dV  = grid->dV[k][j][i];
    rho = d->Vc[RHO][k][j][i];
    prs = d->Vc[PRS][k][j][i];
    vx1 = d->Vc[VX1][k][j][i];
    vx2 = d->Vc[VX2][k][j][i];
    vx3 = d->Vc[VX3][k][j][i];
    
    Temp = (prs/rho) * pow(UNIT_VELOCITY,2.)*(mu*CONST_mp) / CONST_kB; // Kelvin

    tot_vol   += (1-d->Vc[TRC][k][j][i]) * dV; // background wind
    tot_mass  += rho * (1-d->Vc[TRC][k][j][i])* dV;

    v_rms   += (pow(vx1-vx1_avg, 2.) + pow(vx2-vx2_avg, 2.) + pow(vx3-vx3_avg, 2.)) * (1-d->Vc[TRC][k][j][i]) * dV;
    
    if (rho >= sqrt(chi)/3.) {
       cloud_mass += rho * dV; // density based
       cloud_vol  += dV;
    }
    trc_mass += (rho*d->Vc[TRC][k][j][i]) * dV;
    // center of mass of various components
    cmx1_cl += x1[k] * (rho*d->Vc[TRC][k][j][i]) * dV; // along wind
    cmx2_cl += x2[k] * (rho*d->Vc[TRC][k][j][i]) * dV;
    cmx3_cl += x3[k] * (rho*d->Vc[TRC][k][j][i]) * dV;
    vx1_cl  += vx1 * (rho*d->Vc[TRC][k][j][i]) * dV;
    vx2_cl  += vx2 * (rho*d->Vc[TRC][k][j][i]) * dV;
    vx3_cl  += vx3 * (rho*d->Vc[TRC][k][j][i]) * dV;
  }
  /* --- end of main loop ---- */

  /* ------ Parallel data reduction ------ */
  #ifdef PARALLEL
  count = 11;
  tmp = 0;

  sendArray[tmp++] = tot_mass;    sendArray[tmp++] = tot_vol; 
  sendArray[tmp++] = cloud_mass;  sendArray[tmp++] = cloud_vol;   
  sendArray[tmp++] = cmx1_cl;     sendArray[tmp++] = cmx2_cl;     sendArray[tmp++] = cmx3_cl;
  sendArray[tmp++] = vx1_cl;      sendArray[tmp++] = vx2_cl;      sendArray[tmp++] = vx3_cl;
  sendArray[tmp++] = trc_mass;  

  MPI_Allreduce(sendArray, recvArray, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  tmp = 0;
  tot_mass = recvArray[tmp++];    tot_vol = recvArray[tmp++];

  cloud_mass = recvArray[tmp++];  cloud_vol = recvArray[tmp++];
  cmx1_cl = recvArray[tmp++];     cmx2_cl = recvArray[tmp++];     cmx3_cl = recvArray[tmp++];
  vx1_cl = recvArray[tmp++];      vx2_cl = recvArray[tmp++];      vx3_cl = recvArray[tmp++];
  trc_mass = recvArray[tmp++];
  #endif
  /* --- end of parallel data reduction --- */ 

  /* --- Write ascii file "analysis.dat" to disk --- */
  if (prank == 0) {
     cmx1_cl /= cloud_mass;   cmx2_cl /= cloud_mass;    cmx3_cl /= cloud_mass;
     vx1_cl /= cloud_mass;    vx2_cl /= cloud_mass;     vx3_cl /= cloud_mass;
   
     FILE *fp;
     char fname[512];
     static double tpos = -1.0;
     sprintf (fname, "%s/analysis.dat", RuntimeGet()->output_dir);
     if (g_stepNumber == 0) {    /* Open for writing only when we're starting */
       fp = fopen(fname,"w");     /* from beginning */
       fprintf (fp, "# %s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                  "(0) t",        "(1) dt",       "(2) mass",      "(3) TE",
                  "(4) Vx1Avg",   "(5) Vx2Avg",   "(6) Vx3Avg",
                  "(7) Cl_mass",  "(8) Cl_vol",   "(9) CMX1_cl",
                  "(10) CMX2_cl", "(11) CMX3_cl", "(12) Vx1_cl",
                  "(13) Vx2_cl",  "(14) Vx3_cl",  "(15) trc_mass",
                  "(16) leftEdge"
               );
     }
     else {                  /* write if current time t>0*/
       /* Don't duplicate on restart */
       if (tpos < 0.0) {       /* obtain time coordinate of to last written row*/
          char   sline[512];
          fp = fopen(fname,"r");
          while (fgets(sline, 512, fp))  {}
          sscanf(sline, "%lf\n",&tpos);  /* tpos = time of the last written row */
          fclose(fp);
       }
       fp = fopen(fname,"a");
    }
    if (g_time > tpos) {
       fprintf (fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                   g_time,      g_dt,        tot_mass,       tot_TE, 
                   vx1_avg,     vx2_avg,     vx3_avg,
                   cloud_mass,  cloud_vol,   cmx1_cl,
                   cmx2_cl,     cmx3_cl,     vx1_cl,
                   vx2_cl,      vx3_cl,      trc_mass,
                   g_tracerLeftEdge
               );
    }
    fclose(fp);
  }
  /* --- end of writing "analyis.dat" file --- */
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
   int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];
  
  /* wind properties in code units */
  double rhoWind = 1.0; 
  double mach    = g_inputParam[MACH];
  double pWind   = 1.0/(g_gamma*pow(mach, 2.));
  double vWind   = 1.0;  

  /* cloud properties in code units */
  double chi      = g_inputParam[CHI];
  double x_offset = g_inputParam[XOFFSET];
  double Tcl      = g_inputParam[TCL];

  double dummy[4];
  double mu = MeanMolecularWeight((double*)d->Vc, dummy);

  if (side == 0) {    /* -- check solution inside domain -- */
    TOT_LOOP(k,j,i) {
      double temp = (d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i]) * pow(UNIT_VELOCITY,2) * (mu * CONST_mp)/CONST_kB; // Kelvin
      if (temp < Tcl) {
          d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i] * temp) / ( pow(UNIT_VELOCITY,2) * (mu * CONST_mp)/CONST_kB ); 
      }
    }
  }

  if (side == X1_BEG) {  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i) {
        d->Vc[RHO][k][j][i] = rhoWind;
        d->Vc[PRS][k][j][i] = pWind;
        DIM_EXPAND(
        d->Vc[VX1][k][j][i] = vWind-g_vcloud;,
        d->Vc[VX2][k][j][i] = 0.;,
        d->Vc[VX3][k][j][i] = 0.;
        )
      }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){ 
        d->Vc[RHO][k][j][i] = rhoWind;
        d->Vc[PRS][k][j][i] = pWind;
        d->Vc[VX1][k][j][i] = vWind;
        d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][i];
        d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][i];
      }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }


  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
