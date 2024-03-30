/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Module header file for relativistic MHD (RMHD).

  Set label, indexes and basic prototyping for the relativistic 
  MHD module.

  \authors A. Mignone (mignone@to.infn.it)\n
           C. Zanni   (zanni@oato.inaf.it)\n
           G. Mattia
  \date    Dec 03, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */

#if RADIATION
  #include "RMHD/Radiation/radiation.h"
#endif

/* *********************************************************
    Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.  
   ********************************************************* */

#define  RHO 0
#define  MX1 1
#define  MX2 2
#define  MX3 3
#define  BX1 4
#define  BX2 5
#define  BX3 6

#if HAVE_ENERGY
  #define ENG  7
  #define PRS  ENG
#endif

#if DIVB_CONTROL == DIV_CLEANING
    #define PSI_GLM  (ENG + 1)
    #define DIV_COMP 1
#else
  #define DIV_COMP 0
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#if RADIATION
  #define ENR (7 + DIV_COMP + HAVE_ENERGY)
  #define FR1 (ENR + 1)
  #define FR2 (ENR + 2)
  #define FR3 (ENR + 3)
  #define NFLX (11 + DIV_COMP + HAVE_ENERGY)
#else
  #define NFLX (7 + DIV_COMP + HAVE_ENERGY)
#endif

/* *********************************************************
    Label the different waves in increasing order 
    following the number of vector components.

    IMPORTANT: the KPSI_GLMM & KPSI_GLMP modes are 
               present only in the MHD-GLM formulation.
               We keep them at the END of the enumeration
               so we can skip them in unnecessary loops.
               Please do NOT change them !
   ********************************************************* */

enum KWAVES {
 KFASTM, KFASTP, KENTRP

 #if DIVB_CONTROL != DIV_CLEANING
  , KDIVB
 #endif

  , KSLOWM, KSLOWP, KALFVM, KALFVP

 #if DIVB_CONTROL == DIV_CLEANING  
  , KPSI_GLMM, KPSI_GLMP 
 #endif
};

/* ********************************************************************* */
/*! The Map_param structure is used to pass input/output arguments 
    during the conversion from conservative to primitive variables 
    operated by the ConsToPrim() function in the relativistic modules
    (RHD and RMHD).
    The output parameter, rho, W, lor and p, must be set at the end of
    every root-finder: RMHD_[EnergySolve, EntropySolve, PressureFix]().
    Additionally, some of the input parameters must be re-computed in
    EntropySolve() and PressureFix().
   ********************************************************************* */
typedef struct Mao_param_{
 double D;       /**< Lab density       (input). */
 double sigma_c; /**< Conserved entropy (input). */
 double E;       /**< Total energy      (input). */
 double m2;      /**< Square of total momentum (input). */
 double S;       /**< m<dot>B                  (input). */
 double S2;      /**< Square of S              (input). */
 double B2;      /**< Square of magnetic field (input). */

 double rho;     /**< proper density     (output)  */
 double W;       /**< D*h*lor            (output). */
 double lor;     /**< Lorentz factor     (output). */
 double prs;     /**< Thermal pressure   (output). */
} Map_param;

/* ******************************************************
     Vector potential: these labels are and MUST only
     be used in the STARTUP / INIT  functions;
     they're convenient in obtaining a discretization 
     that preserve divB since the beginning.   
   ****************************************************** */

#define   AX1  (NVAR + 1)
#define   AX2  (NVAR + 2)
#define   AX3  (NVAR + 3)

#define AX  AX1  /* backward compatibility */
#define AY  AX2
#define AZ  AX3

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CYLINDRICAL
                               
 #define iVR    VX1
 #define iVZ    VX2
 #define iVPHI  VX3
                               
 #define iMR    MX1
 #define iMZ    MX2
 #define iMPHI  MX3
                               
 #define iBR    BX1
 #define iBZ    BX2
 #define iBPHI  BX3
                               
#endif

#if GEOMETRY == POLAR
                               
 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3
                               
 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3
                               
 #define iBR    BX1
 #define iBPHI  BX2
 #define iBZ    BX3
                               
#endif

#if GEOMETRY == SPHERICAL
                    
 #define iVR     VX1
 #define iVTH    VX2
 #define iVPHI   VX3
                             
 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3
                               
 #define iBR    BX1
 #define iBTH   BX2
 #define iBPHI  BX3

#endif

/* ******************************************************************
    Module-specific symbolic constants (switches)
   ****************************************************************** */

#ifndef RMHD_FAST_EIGENVALUES
  #define RMHD_FAST_EIGENVALUES  NO /**< If set to YES, use approximate (and
                                         faster) expressions when computing the
                                         fast magnetosonic speed, see Sect. 3.3
                                         of Del Zanna,  A&A (2006), 473.
                                         Solutions of quartic equation is avoided
                                         and replace with upper bounds provided by
                                         quadratic equation. */
#endif  

#ifndef RMHD_REDUCED_ENERGY
  #define RMHD_REDUCED_ENERGY    YES  /**< By turning RMHD_REDUCED_ENERGY to YES, 
                                            we let PLUTO evolve the total energy 
                                            minus the mass density contribution. */
#endif 

                                      

/* ---- Function prototyping ----  */

int  ApproximateFastWaves  (double *, double, double, double *);
int  ConsToPrim   (double **, double **, int, int, unsigned char *);
void ConvertTo4vel (double **, int, int);
void ConvertTo3vel (double **, int, int);
void PrimEigenvectors (double *, double, double, double *, double **, double **);
 
void Flux      (const State *, int, int);
void HLL_Speed (const State *, const State *, double *, double *, int, int);
int  MaxSignalSpeed (const State *, double *, double *, int, int);

void PrimToCons   (double **, double **, int, int);
void VelocityLimiter (double *, double *, double *);

int  Magnetosonic (double *vp, double cs2, double h, double *lambda);

Riemann_Solver LF_Solver, HLL_Solver, HLLC_Solver, HLLD_Solver, HLL_Linde_Solver; 
Riemann_Solver GFORCE_Solver, GMUSTA1_Solver;

int  RMHD_EntropySolve (Map_param *);
int  RMHD_EnergySolve  (Map_param *);
int  RMHD_PressureFix  (Map_param *);

#if DIVB_CONTROL == EIGHT_WAVES
 void POWELL_DIVB_SOURCE(const Sweep *, int, int, Grid *);
 void HLL_DIVB_SOURCE (const Sweep *, double **, int, int, Grid *);
#elif DIVB_CONTROL == DIV_CLEANING
 #include "MHD/GLM/glm.h"
#elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
  #include "MHD/CT/ct.h"
#endif

#ifndef NEW_RMHD_FLUXES
  #define NEW_RMHD_FLUXES  NO
#endif
