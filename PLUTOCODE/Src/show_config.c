/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Print useful information about the current computations.

  \author  A. Mignone (mignone@to.infn.it)
  \date    March 20, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
static void CheckConfig();

/* ********************************************************************* */
void ShowConfig (int argc, char *argv[], char *ini_file)
/*!
 *  Write a summary of the selected options
 *
      ___  __   __  ____________ 
     / _ \/ /  / / / /_  __/ __ \
    / ___/ /__/ /_/ / / / / /_/ /
   / /  /____/\____/ /_/  \____/ 
   ============================== 
                                  
 *
 * 
 *                        
 *
 *********************************************************************** */
{
  int  n;
  FILE *fp;
  time_t time_now;
  char  str1[128], str2[128], str3[128], sline[512];
  double *dbl_pnt;  /* For printLoging size of pointer to double */
  int    *int_pnt;  /* For printLoging size of pointer to int */

  CheckConfig();

  printLog ("\n");
  printLog ("   ___  __   __  ____________   \n");
  printLog ("  / _ \\/ /  / / / /_  __/ __ \\ \n");
  printLog (" / ___/ /__/ /_/ / / / / /_/ /  \n");
  printLog ("/_/  /____/\\____/ /_/  \\____/   \n");
  printLog ("=============================    v. %s  \n", PLUTO_VERSION);
  
  printLog ("\n> System:\n\n");

  if ( (fp = fopen("sysconf.out","r")) != NULL){

    while (fscanf (fp, "%s %s %s\n", str1, str2, str3) != EOF) {
      if (!strcmp(str1,"USER")) 
        printLog ("  %s:             %s\n",str1, str3);
      else if (!strcmp(str1,"WORKING_g_dir"))
        printLog ("  %s:      %s\n",str1, str3);
      else if (!strcmp(str1,"SYSTEM_NAME"))
        printLog ("  %s:      %s\n",str1, str3);
      else if (!strcmp(str1,"NODE_NAME"))
        printLog ("  %s:        %s\n",str1, str3);
      else if (!strcmp(str1,"ARCH"))
        printLog ("  %s:             %s\n",str1, str3);
      else if (!strcmp(str1,"BYTE_ORDER"))
        printLog ("  %s:       %s\n\n",str1, str3);
    }
    fclose(fp);

  }else{
    printLog ("  sysconf.out file not found... \n\n");
  }
 
  printLog ("  - Basic data type:\n");
  printLog ("    o sizeof (char)     = %d\n", sizeof(char));
  printLog ("    o sizeof (uchar)    = %d\n", sizeof(unsigned char));
  printLog ("    o sizeof (short)    = %d\n", sizeof(short));
  printLog ("    o sizeof (ushort)   = %d\n", sizeof(unsigned short));
  printLog ("    o sizeof (int)      = %d\n", sizeof(int));
  printLog ("    o sizeof (long)     = %d\n", sizeof(long));
  printLog ("    o sizeof (*int)     = %d\n", sizeof(int_pnt));
  printLog ("    o sizeof (float)    = %d\n", sizeof(float));
  printLog ("    o sizeof (double)   = %d\n", sizeof(double));
  printLog ("    o sizeof (*double)  = %d\n", sizeof(dbl_pnt));
  

  printLog ("\n  - Structure data type:\n");
  printLog ("    o sizeof (cmdLine)    = %d\n", sizeof(cmdLine));
  printLog ("    o sizeof (Data)       = %d\n", sizeof(Data));
  printLog ("    o sizeof (Grid)       = %d\n", sizeof(Grid));
  printLog ("    o sizeof (Float_Vect) = %d\n", sizeof(Float_Vect));
  printLog ("    o sizeof (Image)      = %d\n", sizeof(Image));
  printLog ("    o sizeof (Output)     = %d\n", sizeof(Output));
  printLog ("    o sizeof (RGB)        = %d\n", sizeof(RGB));
  printLog ("    o sizeof (Runtime)    = %d\n", sizeof(Runtime));
  printLog ("    o sizeof (Restart)    = %d\n", sizeof(Restart));
  printLog ("    o sizeof (timeStep)   = %d\n", sizeof(timeStep));
  printLog ("    o sizeof (RBox)       = %d\n", sizeof(RBox));
  printLog ("    o sizeof (State)      = %d\n", sizeof(State));
  printLog ("    o sizeof (Sweep)      = %d\n", sizeof(Sweep));
#if (PARTICLES != NO)
  printLog ("    o sizeof (PARTICLE)   = %d\n", sizeof(Particle));
#endif

  time(&time_now);
  printLog("\n> Local time:       %s\n",asctime(localtime(&time_now)));
      
/* -- printLog command line arguments -- */

  printLog ("> Cmd line args:    ");
  for (n = 1; n < argc; n++) printLog ("%s ",argv[n]);
  printLog ("\n\n");

/* -- printLog problem configuration -- */

  printLog ("> Header configuration:\n\n");

  if (PHYSICS == ADVECTION) printLog ("  PHYSICS:          ADVECTION\n");
  if (PHYSICS == HD)        printLog ("  PHYSICS:          HD\n");
  if (PHYSICS == RHD)       printLog ("  PHYSICS:          RHD\n");
  if (PHYSICS == MHD)       printLog ("  PHYSICS:          MHD    [div.B: ");
  if (PHYSICS == RMHD)      printLog ("  PHYSICS:          RMHD   [div.B: ");
  if (PHYSICS == ResRMHD)   printLog ("  PHYSICS:          ResMHD [div.B: ");
#if (PHYSICS == MHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD)
  #if DIVB_CONTROL == NO
  printLog ("None]\n");
  #elif DIVB_CONTROL == EIGHT_WAVES
    printLog ("Powell's 8wave]\n");
  #elif DIVB_CONTROL == DIV_CLEANING
    #if GLM_EXTENDED == NO 
    printLog ("Divergence Cleaning (GLM)]\n");
    #elif GLM_EXTENDED == YES
    printLog ("Divergence Cleaning (Extended GLM)]\n");
    #endif
    
  #elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
    printLog ("CT/");
    #if CT_EMF_AVERAGE == ARITHMETIC
    printLog ("Ar. average]\n");
    #elif CT_EMF_AVERAGE == CT_CONTACT
    printLog ("CT_CONTACT]\n");
    #elif CT_EMF_AVERAGE == UCT0
    printLog ("UCT0]\n");
    #elif CT_EMF_AVERAGE == UCT_HLL
    printLog ("UCT_HLL]\n");
    #elif CT_EMF_AVERAGE == CT_FLUX
    printLog ("CT_FLUX]\n");
    #elif CT_EMF_AVERAGE == CT_MAXWELL
    printLog ("CT_MAXWELL]\n");
    #elif CT_EMF_AVERAGE == UCT_HLLD
    printLog ("UCT_HLLD]\n");
    #elif CT_EMF_AVERAGE == UCT_GFORCE
    printLog ("UCT_GFORCE]\n");
    #endif
  #endif
#endif

  printLog ("  DIMENSIONS:       %d\n", DIMENSIONS);

  printLog ("  GEOMETRY:         ");
  if (GEOMETRY == CARTESIAN)    printLog ("Cartesian\n");
  if (GEOMETRY == CYLINDRICAL)  printLog ("Cylindrical\n");
  if (GEOMETRY == POLAR)        printLog ("Polar\n");
  if (GEOMETRY == SPHERICAL)    printLog ("Spherical\n");

  printLog ("  BODY_FORCE:       ");
  printLog (BODY_FORCE == NO ? "NO\n":"EXPLICIT\n");

#if COOLING == H2_COOL
  printLog ("  COOLING:          H2_COOL\n");
#elif COOLING == KROME
  printLog ("  COOLING:          KROME\n");
#elif COOLING == MINEq
  printLog ("  COOLING:          MINEq\n");
#elif COOLING == POWER_LAW
  printLog ("  COOLING:          POWER_LAW\n");
#elif COOLING == SNEq
  printLog ("  COOLING:          SNEq\n");
#elif COOLING == TABULATED
  printLog ("  COOLING:          Tabulated\n");
#elif COOLING == TOWNSEND
  printLog ("  COOLING:          Townsend\n");
#elif COOLING == GRACKLE
  printLog ("  COOLING:          Grackle\n");
#endif

  printLog ("  RECONSTRUCTION:   ");
#ifndef FINITE_DIFFERENCE
   if (RECONSTRUCTION == FLAT)        printLog ("Flat");
   if (RECONSTRUCTION == LINEAR)      printLog ("Linear TVD");
   if (RECONSTRUCTION == LimO3)       printLog ("LimO3");
   if (RECONSTRUCTION == WENO3)       printLog ("WENO 3rd order");
   if (RECONSTRUCTION == PARABOLIC)   printLog ("Parabolic");
   if (RECONSTRUCTION == WENOZ)       printLog ("WENOZ (5-th order)");
   if (RECONSTRUCTION == MP5)         printLog ("MP5");
 #ifdef CHAR_LIMITING
   if (CHAR_LIMITING == YES) printLog (" (Characteristic lim)\n");
   else                      printLog (" (Primitive lim)\n");
 #endif
#endif

#ifdef FINITE_DIFFERENCE
  if (RECONSTRUCTION == LIMO3_FD)     printLog ("LimO3 (FD), 3rd order\n");
  if (RECONSTRUCTION == WENO3_FD)     printLog ("WENO3 (FD), 3rd order\n");
  if (RECONSTRUCTION == WENOZ_FD)     printLog ("WENOZ (FD) 5th order\n");
  if (RECONSTRUCTION == MP5_FD)       printLog ("MP5 (FD), 5th order\n");
#endif

  printLog ("  TRACERS:          %d\n", NTRACER);
  printLog ("  VARIABLES:        %d\n", NVAR);
  printLog ("  ENTROPY_SWITCH:   %s\n",(ENTROPY_SWITCH != NO ? "ENABLED":"NO"));
#if PHYSICS == MHD 
  printLog ("  BACKGROUND_FIELD: %s\n",(BACKGROUND_FIELD == YES ? "YES":"NO"));
#endif

  printLog ("  LOADED MODULES:\n");
  #if PHYSICS == MHD
   #ifdef SHEARINGBOX
    printLog ("\n  o [SHEARINGBOX]\n");
    printLog ("     - Order:             %d\n", SB_ORDER);
    printLog ("     - Sym Hydro Flux:    %s\n", 
             (SB_SYMMETRIZE_HYDRO == YES ? "YES":"NO"));
    printLog ("     - Sym Ey:            %s\n", 
             (SB_SYMMETRIZE_EY == YES ? "YES":"NO"));
    printLog ("     - Sym Ez:            %s\n", 
             (SB_SYMMETRIZE_EZ == YES ? "YES":"NO"));
    printLog ("     - Force EMF periods: %s\n", 
             (SB_FORCE_EMF_PERIODS == YES ? "YES":"NO"));
   #endif
  #endif
  #ifdef FARGO
  printLog ("\n  o [FARGO]\n");
  printLog ("     - Order:         %d\n", FARGO_ORDER);
  printLog ("     - Av. Frequency: %d\n", FARGO_NSTEP_AVERAGE);
  #endif
  printLog ("\n");

  printLog ("  ROTATION:         ");
  printLog(ROTATING_FRAME == YES ? "YES\n":"NO\n");

  printLog ("  EOS:              ");
  if      (EOS == IDEAL)        printLog ("Ideal\n");
  else if (EOS == PVTE_LAW)     printLog ("PVTE_LAW\n");
  else if (EOS == BAROTROPIC)   printLog ("Barotropic\n");
  else if (EOS == ISOTHERMAL)   printLog ("Isothermal\n");
  else if (EOS == TAUB)         printLog ("Taub - TM\n");
  else                          printLog ("None\n");

  printLog ("  TIME STEPPING:    ");
  if (TIME_STEPPING == EULER)            printLog ("Euler\n");
  if (TIME_STEPPING == RK2)              printLog ("Runga-Kutta II\n");
  if (TIME_STEPPING == RK3)              printLog ("Runga_Kutta III\n");
  if (TIME_STEPPING == CHARACTERISTIC_TRACING)
                                         printLog ("Characteristic Tracing\n");
#if TIME_STEPPING == HANCOCK
  if (PRIMITIVE_HANCOCK == YES) printLog ("Hancock [Primitive]\n");
  else                          printLog ("Hancock [Conservative]\n");
#endif

  #if PARABOLIC_FLUX != NO
   printLog ("  DIFFUSION TERMS:");
   #if (RESISTIVITY == EXPLICIT) 
    printLog ("  Resistivity  [EXPLICIT]\n");  
   #elif (RESISTIVITY == SUPER_TIME_STEPPING)
    printLog ("  Resistivity  [STS]\n");  
   #elif (RESISTIVITY == RK_LEGENDRE)
    printLog ("  Resistivity  [RKL]\n");
   #endif

   #if (THERMAL_CONDUCTION == EXPLICIT) 
    printLog ("  Thermal Conduction [EXPLICIT]\n");  
   #elif (THERMAL_CONDUCTION == SUPER_TIME_STEPPING)
    printLog ("  Thermal Conduction [STS]\n");  
   #elif (THERMAL_CONDUCTION == RK_LEGENDRE)
    printLog ("  Thermal Conduction [RKL]\n");
   #endif

   #if (VISCOSITY == EXPLICIT) 
    printLog ("  Viscosity [EXPLICIT]\n");  
   #elif (VISCOSITY == SUPER_TIME_STEPPING)
    printLog ("  Viscosity [STS]\n");  
   #elif (VISCOSITY == RK_LEGENDRE)
    printLog ("  Viscosity [RKL]\n");
   #endif
  #endif

  printLog ("\n");

/* -----------------------------------------------------
    Print runtime configuration info (definitions.h 
    and from pluto.ini)
   ----------------------------------------------------- */
/*   
  printLog ("> Header file configuration (definitions.h):\n\n");
  printLog ("  +----------------------------------------------------------\n");
  fp = fopen("definitions.h","r");
  while ( fgets(sline, 512, fp) != NULL ) {
    printLog ("  | %s",sline);
  }
  fclose(fp);
  printLog ("  +---------------------------------------------------------\n\n");
*/
  printLog ("> Runtime configuration (%s):\n\n", ini_file);
  printLog ("  +----------------------------------------------------------\n");
  fp = fopen(ini_file,"r");
  while ( fgets(sline, 512, fp) != NULL ) {
    printLog ("  | %s",sline);
  }
  fclose(fp);
  printLog ("  +---------------------------------------------------------\n");

}

/* ********************************************************************* */
void ShowUnits ()
/*!
 *  Show units when cooling is enabled.
 *
 *
 *********************************************************************** */
{

#if COOLING != NO
  printLog ("> Cooling Module:    ");
  if (COOLING == SNEq)      printLog (" SNEq\n");
  if (COOLING == MINEq)     printLog (" MINEq\n");
  if (COOLING == TABULATED) printLog (" Tabulated\n");
  if (COOLING == TOWNSEND)  printLog (" Townsend\n");
  if (COOLING == H2_COOL)   printLog (" H2_COOL \n");
  if (COOLING == KROME)     printLog (" KROME \n");
  if (COOLING == GRACKLE)   printLog (" GRACKLE \n");
#endif

#if COOLING==GRACKLE
  char grackle_version_info[126];
  grackle_cooling_version_info(grackle_version_info);
  printLog ("> Grackle cooling version: %s\n", grackle_version_info);
#endif

  printLog ("> Normalization Units:\n\n");
  printLog ("  [Density]:      %8.3e (gr/cm^3), %8.3e (1/cm^3)\n",
          UNIT_DENSITY,UNIT_DENSITY/CONST_mp);
  printLog ("  [Pressure]:     %8.3e (dyne/cm^2)\n",
          UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
  printLog ("  [Velocity]:     %8.3e (cm/s)\n",UNIT_VELOCITY);
  printLog ("  [Length]:       %8.3e (cm)\n",UNIT_LENGTH);
  printLog ("  [Temperature]:  %8.3e X (p/rho*mu) (K)\n",KELVIN);
  printLog ("  [Time]:         %8.3e (sec), %8.3e (yrs) \n",
       UNIT_LENGTH/UNIT_VELOCITY, UNIT_LENGTH/UNIT_VELOCITY/86400./365.);
#if PHYSICS == MHD || PHYSICS == RMHD
  printLog ("  [Mag Field]:    %8.3e (Gauss)\n",
           UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY));
#endif

  printLog (" \n");
}

/* ********************************************************************* */
void CheckConfig()
/*
 *
 *
 * Check if the selected configuration is 
 * allowed.
 *
 *
 *********************************************************************** */
{
#if DIMENSIONS == 3 

  #if GEOMETRY  == CYLINDRICAL 
  printLog ("\n! Cylindrical coordinates are only 2D.\n");
  printLog ("! Use polar instead.\n");
  QUIT_PLUTO(1);
  #endif

  #if GEOMETRY == SPHERICAL 
  #if (TIME_STEPPING == HANCOCK) || (TIME_STEPPING == CHARACTERISTIC_TRACING)
  printLog ("\n ! Spherical 3D only works with RK integrators\n");
  QUIT_PLUTO(1);
  #endif
  #endif

#endif
   
}
