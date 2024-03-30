/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Collection of functions for debugging purposes.
  
  \author A. Mignone (mignone@to.infn.it)
  \date   Apr 03, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h" 

/* ********************************************************************* */
int CheckNaN (double **u, int ibeg, int iend, const char *str)
/*!
 * Check whether the array \c u contains Not-a-Number
 *  (NaN). QUIT if true.
 *
 * \param [in]  u      pointer to an array of type u[i][nv]
 * \param [in]  ibeg   starting index 
 * \param [in]  ibeg   ending  index 
 * \param [in]  str    a reference string 
 *
 *********************************************************************** */
{
  int i, nv;

  for (i = ibeg; i <= iend; i++) {
  for (nv = 0; nv < NVAR; nv++) {
    if (isnan(u[i][nv])) {
      printLog ("> CheckNan() [%s]: NaN found", str);
      Where (i,NULL);
      Show (u, i);
      QUIT_PLUTO(1);
    }
  }}
  return 0;
}

/* ********************************************************************* */
int CheckData (Data *data, Grid *grid, const char *str)
/*
 * Check whether the array \c data->Uc or \c data->Vs contain
 * Not-a-Number (NaN), negative density or pressures.
 *
 * \param [in]  data    pointer to PLUTO data strcture
 * \param [in]  grid    pointer to grid structure
 * \param [in]  str    a reference string
 *
 * \return Return 0 if no problem has been found. Otherwise return
 *         1: for nan, 2: for negative density and 3: for negative
 *         pressure or energy.
 *********************************************************************** */
{
  int nv, i,j,k;
  int err = 0;
  char *err_message[] = {" ",
                         "nan found",
                         "negative density",
                         "negative pressure or energy"};

/* --------------------------------------------------------
   1. Check nan for conserved variables
   -------------------------------------------------------- */

  g_dir = IDIR;
  DOM_LOOP(k,j,i){
    g_j = j;
    g_k = k;
    err = 0;
    NVAR_LOOP(nv){
      if (isnan(data->Vc[nv][k][j][i])) err = 1;
    }

    if (data->Vc[RHO][k][j][i] <= 0.0) err = 2;
    #if HAVE_ENERGY
    if (data->Vc[PRS][k][j][i] <= 0.0) err = 3;
    #endif

    if (err > 0){
      printLog ("! CheckData() [%s]: %s\n",str, err_message[err]);
      Where (i,NULL);
      for (nv = 0; nv < NVAR; nv++){
        printLog ("  Vc[%d] = %8.3e\n",nv,data->Vc[nv][k][j][i]);
      }
    }
  }

/* --------------------------------------------------------
   2. Check nan for conserved variables
   -------------------------------------------------------- */

  DOM_LOOP(k,j,i){
    g_j = j;
    g_k = k;
    err = 0;
    NVAR_LOOP(nv){
      if (isnan(data->Uc[k][j][i][nv])) err = 1;
    }

    if (data->Uc[k][j][i][RHO] <= 0.0) err = 2;
    #if HAVE_ENERGY
    if (data->Uc[k][j][i][ENG] <= 0.0) err = 3;
    #endif
    
    if (err){
      printLog ("! CheckData() [%s]: %s\n",str, err_message[err]);
      Where (i,NULL);
      for (nv = 0; nv < NVAR; nv++){
        printLog ("  Uc[%d] = %8.3e\n",nv,data->Uc[k][j][i][nv]);
      }
    }
  }

/* --------------------------------------------------------
   3. Check nan for staggered fields
   -------------------------------------------------------- */

#ifdef STAGGERED_MHD
  TOT_LOOP(k,j,i){
    g_j = j;
    g_k = k;
    err = 0;
    DIM_LOOP(nv){
      if (isnan(data->Vs[nv][k][j][i])) err = 1;
      #if PHYSICS == ResRMHD
      if (isnan(data->Vs[EX1s+nv][k][j][i])) err = 1;
      #endif
    }
    if (err){
      printLog ("! CheckData() [%s]: nan found\n", str);
      Where (i,NULL);
      DIM_LOOP(nv){
        printLog ("! Vs[%d] = %8.3e\n",nv,data->Vs[nv][k][j][i]);
        #if PHYSICS == ResRMHD
        printLog ("! Vs[%d] = %8.3e\n",nv,data->Vs[EX1s+nv][k][j][i]);
        #endif
      }
    }
  }
#endif

  return err;
}

/* ********************************************************************* */
void Show (double **a, int ip)
/*! 
 * Print the component of the array \c a at grid index \c ip  
 *
 *********************************************************************** */
{
  int nv, ix, iy, iz;

  if (g_dir == IDIR) {
    printLog ("X-sweep");
    ix = ip;
    iy = g_j;
    iz = g_k;
  } else if (g_dir == JDIR) {
    printLog ("Y-sweep");
    ix = g_i;
    iy = ip;
    iz = g_k;
  } else if (g_dir == KDIR) {
    printLog ("Z-sweep");
    ix = g_i;
    iy = g_j;
    iz = ip;
  }

  DIM_SELECT( printLog (" (%d)> ", ix);     ,
            printLog (" (%d,%d)> ", ix, iy);  ,
            printLog (" (%d,%d,%d)> ", ix, iy, iz);  )

  for (nv = 0; nv < NVAR; nv++) {
    printLog ("%8.3e  ", a[ip][nv]);
  }
  printLog ("\n");
}

/* ********************************************************************* */
void ShowMatrix(double **A, int n, double eps)
/*!
 * Make a nice printLoging of a 2D square  matrix \c A[0..n-1][0..n-1]
 * Entries with values below eps will display "0.0"
 *
 *********************************************************************** */
{
  int k1,k2;

  printLog ("----------------------------------------------------------------\n");
  for (k1 = 0; k1 < n; k1++){
    for (k2 = 0; k2 < n; k2++){
      printLog ("%12.3e   ", fabs(A[k1][k2]) > eps ? A[k1][k2]:0.0);
    }
    printLog ("\n");
  }
  printLog ("----------------------------------------------------------------\n");
}

/* ********************************************************************* */
void ShowState (double *vc, int primitive)
/*! 
 * Print the components of the array v.
 * Call it as ShowState (vc,1) if vc[] is an array of primitive variables
 * Call it as ShowState (vc,0) if vc[] is an array of conservative variables
 *
 *********************************************************************** */
{
  int nv;

  if (g_dir == IDIR) {
    printLog ("X1-sweep\n");
  } else if (g_dir == JDIR) {
    printLog ("X2-sweep\n");
  } else if (g_dir == KDIR) {
    printLog ("X3-sweep\n");
  }

  if (primitive){
    printLog ("  vc[RHO] = %18.12e;\n", vc[RHO]);
    printLog ("  vc[VX1] = %18.12e;\n", vc[VX1]);
    printLog ("  vc[VX2] = %18.12e;\n", vc[VX2]);
    printLog ("  vc[VX3] = %18.12e;\n", vc[VX3]);
    #if HAVE_ENERGY
    printLog ("  vc[PRS] = %18.12e;\n", vc[PRS]);
    #endif
  
    #if PHYSICS == MHD || PHYSICS == RMHD || PHYSICS == ResRMHD
    printLog ("  vc[BX1] = %18.12e;\n", vc[BX1]);
    printLog ("  vc[BX2] = %18.12e;\n", vc[BX2]);
    printLog ("  vc[BX3] = %18.12e;\n", vc[BX3]);
    #endif
    #if PHYSICS == ResRMHD
    printLog ("  vc[EX1] = %18.12e;\n", vc[EX1]);
    printLog ("  vc[EX2] = %18.12e;\n", vc[EX2]);
    printLog ("  vc[EX3] = %18.12e;\n", vc[EX3]);
    #ifdef CRG
    printLog ("  vc[CRG] = %18.12e;\n", vc[CRG]);
    #endif
    #endif
   
    #ifdef PHI_GLM
    printLog ("  vc[PHI_GLM] = %18.12e;\n", vc[PHI_GLM]);
    #endif
    #ifdef PSI_GLM
    printLog ("  vc[PSI_GLM] = %18.12e;\n", vc[PSI_GLM]);
    #endif
    printLog ("\n");
  }else{
    printLog ("  uc[RHO] = %18.12e;\n", vc[RHO]);
    printLog ("  uc[MX1] = %18.12e;\n", vc[MX1]); 
    printLog ("  uc[MX2] = %18.12e;\n", vc[MX2]); 
    printLog ("  uc[MX3] = %18.12e;\n", vc[MX3]);
    #if HAVE_ENERGY
    printLog ("  uc[ENG] = %18.12e;\n", vc[ENG]);
    #endif
  
    #if PHYSICS == MHD || PHYSICS == RMHD || PHYSICS == ResRMHD
    printLog ("  uc[BX1] = %18.12e;\n", vc[BX1]);
    printLog ("  uc[BX2] = %18.12e;\n", vc[BX2]);
    printLog ("  uc[BX3] = %18.12e;\n", vc[BX3]);
    #endif
    #if PHYSICS == ResRMHD
    printLog ("  uc[EX1] = %18.12e;\n", vc[EX1]);
    printLog ("  uc[EX2] = %18.12e;\n", vc[EX2]);
    printLog ("  uc[EX3] = %18.12e;\n", vc[EX3]);
    #ifdef CRG
    printLog ("  uc[CRG] = %18.12e;\n", vc[CRG]);
    #endif
    #endif
   
    #ifdef PHI_GLM
    printLog ("  uc[PHI_GLM] = %18.12e;\n", vc[PHI_GLM]);
    #endif
    #ifdef PSI_GLM
    printLog ("  uc[PSI_GLM] = %18.12e;\n", vc[PSI_GLM]);
    #endif
    printLog ("\n");
    
  }
}

/* ********************************************************************* */
void ShowVector (double *v, int n)
/*! 
 * Print the first n components of the vector v[]  
 *
 *********************************************************************** */
{
  int k;

  for (k = 0; k < n; k++)  printLog ("%+12.6e  ", v[k]);
  printLog ("\n");
}

/* ********************************************************************* */
void Trace (double xx)
/*!
 * Print a number xx and the number of times it has been called.
 *
 *********************************************************************** */
{
  static int ik;

  printLog ("Trace ------> %f ,  %d\n", xx, ++ik);
}

/* ********************************************************************* */
void Where (int i, Grid *grid)
/*!
 *  Print the location of a particular zone (i,j,k)
 *  in the computational domain.
 *  \note This function must be initialized before using it 
 *        to store grid information. This is done  by calling 
 *        Where(i, grid) the very first time.
 *        Subsequent calls can be then done by simply using 
 *        Where(i,NULL). 
 *
 *********************************************************************** */
{
  int    ii=0, jj=0, kk=0;
  double x1, x2, x3;
  static Grid *grid_copy;

/* --------------------------------------------------
    Keep a local copy of grid for subsequent calls
   -------------------------------------------------- */
 
  if (grid != NULL){
    grid_copy = grid;
    return;
  }

  #ifdef CH_SPACEDIM
   if (g_intStage < 0) return; /* HOT FIX used by CHOMBO
                             (g_intStage = -1) when writing HDF5 file */
  #endif

/* -- ok, proceed normally -- */
  
  if (g_dir == IDIR){
    DIM_EXPAND(ii = i;, jj = g_j;, kk = g_k;)
  }else if (g_dir == JDIR){
    DIM_EXPAND(ii = g_i;, jj = i;, kk = g_k;)
  }else if (g_dir == KDIR){
    DIM_EXPAND(ii = g_i;, jj = g_j;, kk = i;)
  }

  DIM_EXPAND(
    x1 = grid_copy->x[IDIR][ii];  ,
    x2 = grid_copy->x[JDIR][jj];  ,
    x3 = grid_copy->x[KDIR][kk];
  )

  printLog ("  @step = %d (stage = %d);", g_stepNumber, g_intStage);
  DIM_SELECT(
    printLog (" [i = %d], [x1 = %f]", ii, x1);  ,

    printLog (" [i,j = %d, %d], [x1,x2 =  %f, %f]", ii, jj, x1, x2);  ,

    printLog (" [i,j,k = %d, %d, %d], [x1,x2,x3 = %f, %f, %f]", ii, jj, kk,
               x1, x2, x3);
  )

  #ifdef CHOMBO
  printLog (", Level = %d\n", grid_copy->level);
  return;
  #endif
  printLog ("\n");
}
