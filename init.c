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
  g_minCoolingTemp = g_inputParam[TCL];
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

  g_dist_lab = x_offset;

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
  int k, j, i, nv;
  int yr = 365*24*60*60;
  double t_anl;
  static double trc0 = 0.;
  static double trc0_all = 0.;
  static double chi, eta, mach, Tcl, tcc, factor_rho, factor_temp, mu, rho_cl;
  static int first = 0;
  static long int nstep = -1;
  double temperature_cut[] = {1.2, 2.0, 3.0, 5.0, 10.0};
  double rho_cut[] = {1.2, 2.0, 3.0, 5.0, 10.0};

  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  if (first==0) {
    first = 1;
    double oth_mu[4];
    mu   = MeanMolecularWeight((double*)d->Vc, oth_mu);
    chi    = g_inputParam[CHI];
    eta    = g_inputParam[ETA];
    mach   = g_inputParam[MACH];
    Tcl    = g_inputParam[TCL];    

    /* code units */
    tcc      = sqrt(chi);
    factor_rho    = sqrt(chi)/3;
    factor_temp   = sqrt(eta)/3;
    rho_cl   = chi;
  }
  if (g_stepNumber==0){
    double dV;
    DOM_LOOP(k,j,i){
      dV = grid->dV[k][j][i];
      trc0  += d->Vc[RHO][k][j][i]*d->Vc[TRC][k][j][i]*dV;
    }
    #ifdef PARALLEL
    int transfer_size = 1;
    int transfer = 0;
    double sendArray[transfer_size], recvArray[transfer_size];
    sendArray[transfer++] = trc0;
    MPI_Allreduce (sendArray, recvArray, transfer_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    transfer = 0;
    trc0_all = recvArray[transfer++];
    #else
    trc0_all    = trc0;
    #endif
  }
  if (g_stepNumber==0 && trc0_all == 0) {
    printLog("> Analysis(): Check initialization! Likely some error as no cloud tracer has been detected!\n");
    QUIT_PLUTO(1);
  }

  if (trc0_all==0) { // Means we have restarted!
    g_restart = 1;
    FILE *fp;
    char fname[512];
    int dummy;
    sprintf (fname, "%s/restart-analysis.out",RuntimeGet()->output_dir);
    fp = fopen(fname,"r");
    dummy = fscanf(fp, "%lf", &trc0_all);
    dummy = fscanf(fp, "%ld", &nstep);
    dummy = fscanf(fp, "%lf", &t_anl);
    fclose(fp);
  }

  double Tcutoff = (Tcl>1.0e+04)?Tcl:1.0e+04;
  double Tmax    = 1.e8;

  if (g_stepNumber<=nstep && g_time<=(t_anl+0.5*g_anl_dt)) return;
  g_restart = 0;

  double      trc   = 0.,      trc_all    = 0.;
  double mass_dense = 0., mass_dense_all  = 0.;
  double vx_cloud = 0., vy_cloud = 0., vz_cloud = 0.;
  double vx_cloud_all = 0., vy_cloud_all = 0., vz_cloud_all = 0.;

  double mass_cold[(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0]))];
  double mass_cold_all[(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0]))];

  double mass_cloud[(int)(sizeof(rho_cut) / sizeof(rho_cut[0]))];
  double mass_cloud_all[(int)(sizeof(rho_cut) / sizeof(rho_cut[0]))];

  for (i=0; i<(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])); i++){
    mass_cold[i] = 0.;
    mass_cold_all[i] = 0.;
  }

  for (i=0; i<(int)(sizeof(rho_cut) / sizeof(rho_cut[0])); i++){
    mass_cloud[i] = 0.;
    mass_cloud_all[i] = 0.;
  }

  double dV, rByrInj, rho_wind, T_wind, T_gas;
  int cold_indx;
  int cloud_indx;
  rho_wind = 1.0;
  DOM_LOOP(k,j,i){
    dV = grid->dV[k][j][i]; // Cell volume
    trc         += d->Vc[RHO][k][j][i]*d->Vc[TRC][k][j][i]*dV;
    vx_cloud    += d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*d->Vc[TRC][k][j][i]*dV;
    vy_cloud    += d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*d->Vc[TRC][k][j][i]*dV;
    vz_cloud    += d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*d->Vc[TRC][k][j][i]*dV;

    if(d->Vc[RHO][k][j][i] >= (rho_cl/factor_rho))
      mass_dense += d->Vc[RHO][k][j][i]*dV;

    T_wind = MIN(MAX(eta*Tcl, Tcutoff), Tmax);

    T_gas = (d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*pow(UNIT_VELOCITY,2)*(CONST_mp*mu)/CONST_kB;
    for (cloud_indx=0; cloud_indx<(int)(sizeof(rho_cut) / sizeof(rho_cut[0])); cloud_indx++){
        if (d->Vc[RHO][k][j][i] >= (rho_wind*rho_cut[cloud_indx])){
          if( T_gas <= (factor_temp*Tcl) )
            mass_cloud[cloud_indx] += d->Vc[RHO][k][j][i]*dV;
        }
    }
    if( T_gas <= (factor_temp*Tcl) ){ 
      for (cold_indx=0; cold_indx<(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])); cold_indx++){
          if( T_gas <= (T_wind/temperature_cut[cold_indx]) )
            mass_cold[cold_indx] += d->Vc[RHO][k][j][i]*dV;
      }
      /* if (d->Vc[RHO][k][j][i]>(rho_cl/sqrt(chi))) mass_cold += d->Vc[RHO][k][j][i]*dV; */
      /* if (d->Vc[RHO][k][j][i] >= (rho_cl/factor_rho)) mass_cold += d->Vc[RHO][k][j][i]*dV; */
    }
  }

  #ifdef PARALLEL
  int transfer_size = 5 + (int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])) + (int)(sizeof(rho_cut) / sizeof(rho_cut[0]));
  int transfer = 0;
  double sendArray[transfer_size], recvArray[transfer_size];
  sendArray[transfer++] = trc; sendArray[transfer++] = mass_dense;
  for (cold_indx=0; cold_indx<(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])); cold_indx++) {
    sendArray[transfer++] = mass_cold[cold_indx];
  }
  for (cloud_indx=0; cloud_indx<(int)(sizeof(rho_cut) / sizeof(rho_cut[0])); cloud_indx++) {
    sendArray[transfer++] = mass_cloud[cloud_indx];
  }
  sendArray[transfer++] = vx_cloud; sendArray[transfer++] = vy_cloud; sendArray[transfer++] = vz_cloud;
  MPI_Allreduce (sendArray, recvArray, transfer_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // TODO: Replace this with Allreduce to improve on communication overhead
  transfer = 0;
  trc_all = recvArray[transfer++]; mass_dense_all = recvArray[transfer++];
  for (cold_indx=0; cold_indx<(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])); cold_indx++) {
    mass_cold_all[cold_indx] = recvArray[transfer++];
  }
  for (cloud_indx=0; cloud_indx<(int)(sizeof(rho_cut) / sizeof(rho_cut[0])); cloud_indx++) {
    mass_cloud_all[cloud_indx] = recvArray[transfer++];
  }
  vz_cloud_all = recvArray[transfer++]; vy_cloud_all = recvArray[transfer++]; vz_cloud_all = recvArray[transfer++];

  #else
  trc_all    = trc;
  mass_dense_all = mass_dense;

  for (cold_indx=0; cold_indx<(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])); cold_indx++) {
    mass_cold_all[cold_indx] = mass_cold[cold_indx];
  }
  for (cloud_indx=0; cloud_indx<(int)(sizeof(rho_cut) / sizeof(rho_cut[0])); cloud_indx++) {
    mass_cloud_all[cloud_indx] = mass_cloud[cloud_indx];
  }
  vx_cloud_all    = vx_cloud;
  vy_cloud_all    = vy_cloud;
  vz_cloud_all    = vz_cloud;
  #endif
  vx_cloud_all = vx_cloud_all/trc_all;
  vy_cloud_all = vy_cloud_all/trc_all;
  vz_cloud_all = vz_cloud_all/trc_all;
  trc_all     = trc_all/trc0_all; // trc0_all is M_cloud, ini
  mass_dense_all = mass_dense_all/trc0_all;
  for (cold_indx=0; cold_indx<(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])); cold_indx++) {
    mass_cold_all[cold_indx] = mass_cold_all[cold_indx]/trc0_all;
  }
  for (cloud_indx=0; cloud_indx<(int)(sizeof(rho_cut) / sizeof(rho_cut[0])); cloud_indx++) {
    mass_cloud_all[cloud_indx] = mass_cloud_all[cloud_indx]/trc0_all;
  }

  double v_cloud = sqrt(vx_cloud_all*vx_cloud_all + vy_cloud_all*vy_cloud_all + vz_cloud_all*vz_cloud_all);
  g_dist_lab += (vx_cloud_all+g_vcloud)*g_anl_dt;

  /* ---- Write ascii file "analysis.dat" to disk ---- */
  if (prank == 0){
    char fname[512];
    char buffer1[128], buffer2[128];
    sprintf(buffer1, "M(rho>rho_cl/%.1f)/M0", factor_rho);
    sprintf(buffer2, "M(T<%.1f*T_cl)/M0", factor_temp);
    char *cold_header[(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0]))];
    char *cloud_header[(int)(sizeof(rho_cut) / sizeof(rho_cut[0]))];

    for (cold_indx=0; cold_indx<(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])); cold_indx++) {
      char *dummy1 = (char *)malloc(256*sizeof(char));
      char *dummy2 = (char *)malloc(256*sizeof(char));
      cold_header[cold_indx] = (char *)malloc(256*sizeof(char));
      strcpy(dummy1, buffer2);
      sprintf(dummy2, " [T<T_w/%.1f]", temperature_cut[cold_indx]);
      strcat(dummy1, dummy2);
      strcpy(cold_header[cold_indx], dummy1);
    }

    for (cloud_indx=0; cloud_indx<(int)(sizeof(rho_cut) / sizeof(rho_cut[0])); cloud_indx++) {
      cloud_header[cloud_indx] = (char *)malloc(256*sizeof(char));
      sprintf(cloud_header[cloud_indx], "M (rho>=%.1f rho_w)/M0", rho_cut[cloud_indx]);
    }

    FILE *fp;
    sprintf (fname, "%s/analysis.dat",RuntimeGet()->output_dir);
    if (g_stepNumber == 0){ /* Open for writing only when we are starting */
      fp = fopen(fname,"w"); /* from beginning */
      fprintf (fp,"# %s\t=\t%.5e\n", "tcc (code)", tcc);
      fprintf (fp,"# %s\t=\t%.5e\n", "vwind_asymp (code)", UNIT_VELOCITY );
      // Header
      fprintf (fp,"# (1)%s\t\t(2)%s\t(3)%s\t\t(4)%s\t\t(5)%s\t\t",
               "time (code)", "g_dist_lab (code)", "v_cloud (code)", "trc/trc0", buffer1); //, buffer2, "dt (code)");
      int cont = 5;
      for (cold_indx=0; cold_indx<(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])); cold_indx++) {
        fprintf(fp,"(%d)%s\t\t", ++cont, cold_header[cold_indx]);
      }
      for (cloud_indx=0; cloud_indx<(int)(sizeof(rho_cut) / sizeof(rho_cut[0])); cloud_indx++) {
        fprintf(fp,"(%d)%s\t\t", ++cont, cloud_header[cloud_indx]);
      }
      fprintf (fp, "(%d)dt (code)\n", ++cont);
      fclose(fp);
    }
    /* Append numeric data */
    fp = fopen(fname,"a");
    fprintf (fp, "%12.6e\t\t%12.6e\t\t%12.6e\t\t%12.6e\t\t%12.6e\t\t\t",
             g_time, g_dist_lab, v_cloud, trc_all, mass_dense_all);
    for (cold_indx=0; cold_indx<(int)(sizeof(temperature_cut) / sizeof(temperature_cut[0])); cold_indx++) {
      fprintf (fp, "%12.6e\t\t\t", mass_cold_all[cold_indx]);
    }
    for (cloud_indx=0; cloud_indx<(int)(sizeof(rho_cut) / sizeof(rho_cut[0])); cloud_indx++) {
      fprintf (fp, "%12.6e\t\t\t", mass_cloud_all[cloud_indx]);
    }
    fprintf (fp, "%12.6e\n", g_dt);
    fclose(fp);

    /* Write restart file */
    //printLog("Step %d: Writing Analysis restart!\n", g_stepNumber);
    FILE *frestart;
    sprintf (fname, "%s/restart-analysis.out", RuntimeGet()->output_dir);
    frestart = fopen(fname,"w");
    fprintf (frestart,"%lf\n", trc0_all);
    fprintf(frestart,"%ld\n", g_stepNumber);
    fprintf(frestart,"%lf\n", g_time);
    fclose(frestart);
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
