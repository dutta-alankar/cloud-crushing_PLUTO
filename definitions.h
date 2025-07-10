#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        GRACKLE
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            8

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO
#define  INTERNAL_BOUNDARY              YES
#define  SHOW_TIMING                    NO
#define  SHOW_TIME_STEPS                YES
#define  BOOST                          YES 

/* -- user-defined parameters (labels) -- */

#define  CHI                            0
#define  ETA                            1
#define  MACH                           2
#define  TCOOL_TCC                      3
#define  TCL                            4
#define  XOFFSET                        5
#define  ZMET_CL                        6
#define  ZMET_W                         7

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   1.0472496217017487e-27
#define  UNIT_LENGTH                    1.2804349697480817e+20
#define  UNIT_VELOCITY                  22234757.8993876

/* [End] user-defined constants (do not change this line) */
#define  MULTIPLE_LOG_FILES             YES
#define  VERBOSE                        NO
