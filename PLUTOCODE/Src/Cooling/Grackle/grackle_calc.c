/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute right hand side for Grackle cooling

  \authors A. Dutta (alankard@mpa-garching.mpg.de)\n

 \b References

  \date   May 13, 2025
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "grackle.h"

void grackle_cooling_version_info (char *version) {
    grackle_version gversion = get_grackle_version();
    strcpy(version, gversion.version);
}

void finalize_grackle () {
    free_chemistry_data();
}

void normalize_ions_grackle (const Data *d, const chemistry_data *grackle_config_data, int i, int j, int k) {
    // normalize the ion fractions at cell center
    double norm_H = 0., norm_He = 0.;
    if (grackle_config_data->primordial_chemistry >= 1) {
        norm_H  = d->Vc[X_HI][k][j][i] + d->Vc[X_HII][k][j][i];
        norm_He = d->Vc[Y_HeI][k][j][i] + d->Vc[Y_HeII][k][j][i] + d->Vc[Y_HeIII][k][j][i];
    }
    if (grackle_config_data->primordial_chemistry >= 2) {
        norm_H  += (d->Vc[X_HM][k][j][i] + d->Vc[X_H2I][k][j][i] + d->Vc[X_H2II][k][j][i]);
    }
    if (grackle_config_data->primordial_chemistry >= 3) {
        norm_H  += (d->Vc[X_DI][k][j][i] + d->Vc[X_DII][k][j][i] + d->Vc[X_HDI][k][j][i]);
    }
    // norm_H  /= grackle_config_data->HydrogenFractionByMass;
    // norm_He /= (1-grackle_config_data->HydrogenFractionByMass);
    if (grackle_config_data->primordial_chemistry >= 1) {
        d->Vc[X_HI][k][j][i]    /= norm_H;
        d->Vc[X_HII][k][j][i]   /= norm_H;
        d->Vc[Y_HeI][k][j][i]   /= norm_He;
        d->Vc[Y_HeII][k][j][i]  /= norm_He;
        d->Vc[Y_HeIII][k][j][i] /= norm_He;
        d->Vc[X_HM][k][j][i]    /= norm_H;
        d->Vc[X_H2I][k][j][i]   /= norm_H;
        d->Vc[X_H2II][k][j][i]  /= norm_H;
        d->Vc[X_DI][k][j][i]    /= norm_H;
        d->Vc[X_DII][k][j][i]   /= norm_H;
        d->Vc[X_HDI][k][j][i]   /= norm_H;
    }
}

void call_grackle_equil (const Data *d, Grid *grid) {
    double time = 13.6*1.0e+09*365*24*60*60/(UNIT_LENGTH/UNIT_VELOCITY); // Age of universe
    call_grackle(d, time, NULL, grid, 0, 0, 0, 0);
}

void call_grackle_equil_by_cell (const Data *d, Grid *grid, int i, int j, int k) {
    double time = 13.6*1.0e+09*365*24*60*60/(UNIT_LENGTH/UNIT_VELOCITY); // Age of universe
    call_grackle(d, time, NULL, grid, 1, i, j, k);
}

void call_grackle (const Data *d, double dt, timeStep *Dts, Grid *grid, int one_cell, int cell_i, int cell_j, int cell_k)
/*!
 * Integrate cooling and reaction source terms.
 *
 * \param [in,out]  d     pointer to Data structure
 * \param [in]     dt     the time step to be taken
 * \param [out]    Dts    pointer to the Time_Step structure (do equilibrium if called with NULL)
 * \param [in]     grid   pointer to an array of Grid structures
 * \param[in]      one_cell 1 means do only one cell
 *
 *********************************************************************** */
{
    int i, j, k, id;
    static int once = 0;
    static grackle_field_data grackle_chemistry_fields;
    static chemistry_data *grackle_config_data;
    static timeStep *Dts_prev_call = NULL;
    static int once_cell_prev = 1;

    // This ensures that everything in Grackle is refreshed if one cell or equilibrium was called earlier
    if (Dts_prev_call==NULL || once_cell_prev==1) once = 0;
    Dts_prev_call = Dts;
    once_cell_prev = one_cell;
    
    if (one_cell==1) {
        id = 0;
        once = 0;
    }
    if (once==0 && grackle_chemistry_fields.density!=NULL) {
        free_chemistry_data();
    }
    // if (once==0 && grackle_chemistry_fields.density!=NULL) free_chemistry_data();
    /*********************************************************************
    / Initial setup of units and chemistry objects.
    / This should be done at simulation start.
    *********************************************************************/

    // Set initial redshift (for internal units).
    double initial_redshift = 0.;

    if (once==0) {
        // Check the consistency
        if (gr_check_consistency() != GR_SUCCESS) {
          printLog("call_grackle(): Error in gr_check_consistency.\n");
          QUIT_PLUTO(1);
        }
    }
    // Enable output
    grackle_verbose = g_grackle_params.grackle_verbose;

    // First, set up the units system.
    // These are conversions from code units to cgs.
    code_units grackle_code_units;
    grackle_code_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
    grackle_code_units.density_units = UNIT_DENSITY;
    grackle_code_units.length_units = UNIT_LENGTH;
    grackle_code_units.time_units = UNIT_LENGTH/UNIT_VELOCITY;
    grackle_code_units.a_units = 1.0; // units for the expansion factor
    // Set expansion factor to 1 for non-cosmological simulation.
    grackle_code_units.a_value = 1. / (1. + initial_redshift) / grackle_code_units.a_units;
    set_velocity_units(&grackle_code_units);

    // Second, create a chemistry object for parameters.  This needs to be a pointer.
    if (once==0) {
        if (grackle_config_data==NULL) grackle_config_data = malloc(sizeof(chemistry_data));
        if (set_default_chemistry_parameters(grackle_config_data) == 0) {
            printLog("call_grackle(): Error in set_default_chemistry_parameters.\n");
            QUIT_PLUTO(1);
        }
    }
    // Set parameter values for chemistry.
    // Access the parameter storage with the struct you've created
    // or with the grackle_config_data pointer declared in grackle.h (see further below).
    grackle_config_data->max_iterations = 100000000;
    grackle_config_data->Gamma = g_gamma;
    grackle_config_data->use_grackle = 1;            // chemistry on
    if (Dts!=NULL)
        grackle_config_data->with_radiative_cooling = 1; // cooling on
    else
        grackle_config_data->with_radiative_cooling = 0; // cooling off --> equilibrium 
    grackle_config_data->primordial_chemistry = g_grackle_params.grackle_primordial_chemistry;   // molecular network with H, He, D
    grackle_config_data->dust_chemistry = g_grackle_params.grackle_dust_chemistry;         // dust processes
    grackle_config_data->metal_cooling = g_grackle_params.grackle_metal_cooling;          // metal cooling on
    grackle_config_data->UVbackground = g_grackle_params.grackle_UVbackground;           // UV background on
    grackle_config_data->grackle_data_file = g_grackle_params.grackle_data_file; // data file
    grackle_config_data->use_temperature_floor = g_grackle_params.grackle_use_temperature_floor;  // switch on a scalar temperature floor
    if (g_grackle_params.grackle_temperature_floor_scalar>0)
        grackle_config_data->temperature_floor_scalar = (g_grackle_params.grackle_temperature_floor_scalar>0)?g_grackle_params.grackle_temperature_floor_scalar:g_minCoolingTemp;  // temperature floor

    // Finally, initialize the chemistry object.
    if (once==0) {
        if (initialize_chemistry_data(&grackle_code_units) == 0) {
            printLog("call_grackle(): Error in initialize_chemistry_data.\n");
            QUIT_PLUTO(1);
        }
    }
    double tiny_number = 1.e-20;
    // Create struct for storing grackle field data
    if (once==0) {
        if (grackle_chemistry_fields.grid_dimension!=NULL) FreeArray1D(grackle_chemistry_fields.grid_dimension);
        if (grackle_chemistry_fields.grid_start!=NULL) FreeArray1D(grackle_chemistry_fields.grid_start);
        if (grackle_chemistry_fields.grid_end!=NULL) FreeArray1D(grackle_chemistry_fields.grid_end);
        if (grackle_chemistry_fields.density!=NULL) FreeArray1D(grackle_chemistry_fields.density);
        if (grackle_chemistry_fields.internal_energy!=NULL) FreeArray1D(grackle_chemistry_fields.internal_energy);
        if (grackle_chemistry_fields.x_velocity!=NULL) FreeArray1D(grackle_chemistry_fields.x_velocity);
        if (grackle_chemistry_fields.y_velocity!=NULL) FreeArray1D(grackle_chemistry_fields.y_velocity);
        if (grackle_chemistry_fields.z_velocity!=NULL) FreeArray1D(grackle_chemistry_fields.z_velocity);
        if (grackle_config_data->primordial_chemistry >= 1) {
            if (grackle_chemistry_fields.HI_density!=NULL) FreeArray1D(grackle_chemistry_fields.HI_density);
            if (grackle_chemistry_fields.HII_density!=NULL) FreeArray1D(grackle_chemistry_fields.HII_density);
            if (grackle_chemistry_fields.HeI_density!=NULL) FreeArray1D(grackle_chemistry_fields.HeI_density);
            if (grackle_chemistry_fields.HeII_density!=NULL) FreeArray1D(grackle_chemistry_fields.HeII_density);
            if (grackle_chemistry_fields.HeIII_density!=NULL) FreeArray1D(grackle_chemistry_fields.HeIII_density);
            if (grackle_chemistry_fields.e_density!=NULL) FreeArray1D(grackle_chemistry_fields.e_density);
	}
        if (grackle_config_data->primordial_chemistry >= 2) {
            if (grackle_chemistry_fields.HM_density!=NULL) FreeArray1D(grackle_chemistry_fields.HM_density);
            if (grackle_chemistry_fields.H2I_density!=NULL) FreeArray1D(grackle_chemistry_fields.H2I_density);
            if (grackle_chemistry_fields.H2II_density!=NULL) FreeArray1D(grackle_chemistry_fields.H2II_density);
	}
        if (grackle_config_data->primordial_chemistry >= 3) {
            if (grackle_chemistry_fields.DI_density!=NULL) FreeArray1D(grackle_chemistry_fields.DI_density);
            if (grackle_chemistry_fields.DII_density!=NULL) FreeArray1D(grackle_chemistry_fields.DII_density);
            if (grackle_chemistry_fields.HDI_density!=NULL) FreeArray1D(grackle_chemistry_fields.HDI_density);
	}
        if (grackle_config_data->metal_cooling == 1) 
            if (grackle_chemistry_fields.metal_density!=NULL) FreeArray1D(grackle_chemistry_fields.metal_density);
        if (grackle_chemistry_fields.volumetric_heating_rate!=NULL) FreeArray1D(grackle_chemistry_fields.volumetric_heating_rate);
        if (grackle_chemistry_fields.specific_heating_rate!=NULL) FreeArray1D(grackle_chemistry_fields.specific_heating_rate);
        if (grackle_chemistry_fields.RT_HI_ionization_rate!=NULL) FreeArray1D(grackle_chemistry_fields.RT_HI_ionization_rate);
        if (grackle_chemistry_fields.RT_HeI_ionization_rate!=NULL) FreeArray1D(grackle_chemistry_fields.RT_HeI_ionization_rate);
        if (grackle_chemistry_fields.RT_HeII_ionization_rate!=NULL) FreeArray1D(grackle_chemistry_fields.RT_HeII_ionization_rate);
        if (grackle_chemistry_fields.RT_H2_dissociation_rate!=NULL) FreeArray1D(grackle_chemistry_fields.RT_H2_dissociation_rate);
        if (grackle_chemistry_fields.RT_heating_rate!=NULL) FreeArray1D(grackle_chemistry_fields.RT_heating_rate);
        gr_initialize_field_data(&grackle_chemistry_fields);
    }

    if (once==0) {
        // Set grid dimension and size.
        // grid_start and grid_end are used to ignore ghost zones.
        // int field_size = 1;
        grackle_chemistry_fields.grid_rank = 3;
        if (grackle_chemistry_fields.grid_dimension==NULL)
            grackle_chemistry_fields.grid_dimension = ARRAY_1D(grackle_chemistry_fields.grid_rank, int);
        if (grackle_chemistry_fields.grid_start==NULL)    
            grackle_chemistry_fields.grid_start = ARRAY_1D(grackle_chemistry_fields.grid_rank, int);
        if (grackle_chemistry_fields.grid_end==NULL)    
            grackle_chemistry_fields.grid_end = ARRAY_1D(grackle_chemistry_fields.grid_rank, int);
        grackle_chemistry_fields.grid_dx = -1; // MAX(MAX(grid->dx[KDIR][1], grid->dx[JDIR][1]), grid->dx[IDIR][1]); // used only for H2 self-shielding approximation
    }
    
    for (i = 0; i < 3; i++) {
        grackle_chemistry_fields.grid_dimension[i] = (one_cell==0)?grid->np_int[i]:1; // all cells not including ghost zones (local).
        grackle_chemistry_fields.grid_start[i] = 0;
        grackle_chemistry_fields.grid_end[i] = (one_cell==0)?(grid->np_int[i]-1):0;
    }
    // my_fields.grid_dimension[0] = field_size;
    // my_fields.grid_end[0] = field_size - 1;

    if (once==0) {
        grackle_chemistry_fields.density         = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        grackle_chemistry_fields.internal_energy = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        grackle_chemistry_fields.x_velocity      = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        grackle_chemistry_fields.y_velocity      = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        grackle_chemistry_fields.z_velocity      = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        if (grackle_config_data->primordial_chemistry >= 1) {
            grackle_chemistry_fields.HI_density      = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
            grackle_chemistry_fields.HII_density     = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
            grackle_chemistry_fields.HeI_density     = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
            grackle_chemistry_fields.HeII_density    = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
            grackle_chemistry_fields.HeIII_density   = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
            grackle_chemistry_fields.e_density       = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
	}
        if (grackle_config_data->primordial_chemistry >= 2) {
            grackle_chemistry_fields.HM_density      = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
            grackle_chemistry_fields.H2I_density     = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
            grackle_chemistry_fields.H2II_density    = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
	}
        if (grackle_config_data->primordial_chemistry >= 3) {
            grackle_chemistry_fields.DI_density      = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
            grackle_chemistry_fields.DII_density     = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
            grackle_chemistry_fields.HDI_density     = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
	}
        if (grackle_config_data->metal_cooling == 1)
            grackle_chemistry_fields.metal_density   = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        // volumetric heating rate (provide in units [erg s^-1 cm^-3])
        grackle_chemistry_fields.volumetric_heating_rate = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        // specific heating rate (provide in units [egs s^-1 g^-1]
        grackle_chemistry_fields.specific_heating_rate = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);

        // radiative transfer ionization / dissociation rate fields (provide in units [1/s])
        grackle_chemistry_fields.RT_HI_ionization_rate = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        grackle_chemistry_fields.RT_HeI_ionization_rate = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        grackle_chemistry_fields.RT_HeII_ionization_rate = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        grackle_chemistry_fields.RT_H2_dissociation_rate = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
        // radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
        grackle_chemistry_fields.RT_heating_rate = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
    }
    
    // set temperature units
    double temperature_units = get_temperature_units(&grackle_code_units);

    int counter = 0;
    DOM_LOOP (k, j, i) {
        if (one_cell==0)
            id = (k-grid->lbeg[KDIR]) * grid->np_int[JDIR] * grid->np_int[IDIR] + (j-grid->lbeg[JDIR]) * grid->np_int[IDIR] + (i-grid->lbeg[IDIR]);
        if (one_cell==1) {
            if (counter>0) continue;
            i = cell_i;
            j = cell_j;
            k = cell_k;
        }
        // printLog("DEBUG: (k, j, i, id) = (%d, %d, %d, %d) \n", k, j, i, id);
        if (once==0 || one_cell==1) normalize_ions_grackle(d, grackle_config_data, i, j, k);
	    grackle_chemistry_fields.density[id] = (gr_float)d->Vc[RHO][k][j][i];
        if (grackle_config_data->primordial_chemistry >= 1) {
            grackle_chemistry_fields.HI_density[id] = (gr_float)(d->Vc[X_HI][k][j][i]) * grackle_config_data->HydrogenFractionByMass * grackle_chemistry_fields.density[id];
            grackle_chemistry_fields.HII_density[id] = (gr_float)(d->Vc[X_HII][k][j][i]) * grackle_config_data->HydrogenFractionByMass * grackle_chemistry_fields.density[id];
            grackle_chemistry_fields.HeI_density[id] = (gr_float)(d->Vc[Y_HeI][k][j][i]) * (1-grackle_config_data->HydrogenFractionByMass)  * grackle_chemistry_fields.density[id];
            grackle_chemistry_fields.HeII_density[id] = (gr_float)(d->Vc[Y_HeII][k][j][i]) * (1-grackle_config_data->HydrogenFractionByMass) * grackle_chemistry_fields.density[id];
            grackle_chemistry_fields.HeIII_density[id] = (gr_float)(d->Vc[Y_HeIII][k][j][i]) * (1-grackle_config_data->HydrogenFractionByMass) * grackle_chemistry_fields.density[id];
            grackle_chemistry_fields.e_density[id] = (grackle_chemistry_fields.HII_density[id] + (grackle_chemistry_fields.HeII_density[id] + 2*grackle_chemistry_fields.HeIII_density[id])/4);  
	        // grackle_chemistry_fields.e_density[id] = (gr_float)d->Vc[elec][k][j][i];
	    // normalization: see https://grackle.readthedocs.io/en/latest/Interaction.html#density-note
        }
        if (grackle_config_data->primordial_chemistry >= 2) {
            grackle_chemistry_fields.HM_density[id] = (gr_float)(d->Vc[X_HM][k][j][i]) * grackle_config_data->HydrogenFractionByMass * grackle_chemistry_fields.density[id];
            grackle_chemistry_fields.H2I_density[id] = (gr_float)(d->Vc[X_H2I][k][j][i]) * grackle_config_data->HydrogenFractionByMass * grackle_chemistry_fields.density[id];
            grackle_chemistry_fields.H2II_density[id] = (gr_float)(d->Vc[X_H2II][k][j][i]) * grackle_config_data->HydrogenFractionByMass * grackle_chemistry_fields.density[id];
        }
        if (grackle_config_data->primordial_chemistry >= 3) {
            grackle_chemistry_fields.DI_density[id] = (gr_float)(d->Vc[X_DI][k][j][i]) * grackle_config_data->HydrogenFractionByMass * grackle_chemistry_fields.density[id];
            grackle_chemistry_fields.DII_density[id] = (gr_float)(d->Vc[X_DII][k][j][i]) * grackle_config_data->HydrogenFractionByMass * grackle_chemistry_fields.density[id];
            grackle_chemistry_fields.HDI_density[id] = (gr_float)(d->Vc[X_HDI][k][j][i]) * grackle_config_data->HydrogenFractionByMass * grackle_chemistry_fields.density[id];
        }
        // solar metallicity
        if (grackle_config_data->metal_cooling == 1) 
            grackle_chemistry_fields.metal_density[id] = (gr_float)(d->Vc[Z_MET][k][j][i]) * grackle_config_data->SolarMetalFractionByMass * grackle_chemistry_fields.density[id];

        grackle_chemistry_fields.x_velocity[id] = (gr_float)d->Vc[VX1][k][j][i];
        grackle_chemistry_fields.y_velocity[id] = (gr_float)d->Vc[VX2][k][j][i];
        grackle_chemistry_fields.z_velocity[id] = (gr_float)d->Vc[VX3][k][j][i];

        // initilize specific internal thermal energy
        grackle_chemistry_fields.internal_energy[id] = (gr_float)((d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])/(g_gamma-1));
        /*
        if (grackle_config_data->with_radiative_cooling != 0)
            grackle_chemistry_fields.internal_energy[id] = (gr_float)((d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])/(g_gamma-1));
        else
            grackle_chemistry_fields.internal_energy[id] = (gr_float)(d->Vgrac[TEMP][k][j][i]/(temperature_units*(g_gamma-1)));
        */
        // printLog("DEBUG: (k,j,i)=(%d,%d,%d) T = %e, P/rho = %e  P=%e  rho=%e\n", k, j, i, d->Vgrac[TEMP][k][j][i], (gr_float)((d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])/(g_gamma-1)), d->Vc[PRS][k][j][i], d->Vc[RHO][k][j][i] );

        grackle_chemistry_fields.volumetric_heating_rate[id] = 0.0;
        grackle_chemistry_fields.specific_heating_rate[id] = 0.0;

        grackle_chemistry_fields.RT_HI_ionization_rate[id] = 0.0;
        grackle_chemistry_fields.RT_HeI_ionization_rate[id] = 0.0;
        grackle_chemistry_fields.RT_HeII_ionization_rate[id] = 0.0;
        grackle_chemistry_fields.RT_H2_dissociation_rate[id] = 0.0;
        grackle_chemistry_fields.RT_heating_rate[id] = 0.0;
        counter++;
    }

    /*********************************************************************
    / Calling the chemistry solver
    / These routines can now be called during the simulation.
    *********************************************************************/

    // printLog("Evolving chemistry...\n");
    if (solve_chemistry(&grackle_code_units, &grackle_chemistry_fields, dt) == 0) {
        printLog("call_grackle(): Error in solve_chemistry.\n");
        QUIT_PLUTO(1);
    }

    // Calculate cooling time.
    static gr_float *cooling_time;
    if (once==0 && Dts!=NULL) {
        if (cooling_time!=NULL) FreeArray1D(cooling_time);
        cooling_time = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
    }
    if (Dts!=NULL) {
        if (calculate_cooling_time(&grackle_code_units, &grackle_chemistry_fields,
                                cooling_time) == 0) {
            printLog("call_grackle(): Error in calculate_cooling_time.\n");
            QUIT_PLUTO(1);
        }
        double cool_time_min = 1.0e+30;
        DOM_LOOP(k, j, i) {
            id = (k-grid->lbeg[KDIR]) * grid->np_int[JDIR] * grid->np_int[IDIR] + (j-grid->lbeg[JDIR]) * grid->np_int[IDIR] + (i-grid->lbeg[IDIR]);
            cool_time_min = (fabs(cooling_time[id])<cool_time_min)?fabs(cooling_time[id]):cool_time_min;
        }
        Dts->dt_cool = cool_time_min;
        // printLog("> step %d: cooling_time = %24.16g Myr (%e code)\n", g_stepNumber, cool_time_min * grackle_code_units.time_units/ (1.0e+06*365*24*60*60), cool_time_min );
    }
    // Calculate temperature in K.
    static gr_float *temperature;
    if (once==0) {
        if (temperature!=NULL) FreeArray1D(temperature);
        temperature = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
    }
    if (calculate_temperature(&grackle_code_units, &grackle_chemistry_fields,
                              temperature) == 0) {
        printLog("call_grackle(): Error in calculate_temperature.\n");
        QUIT_PLUTO(1);
    }

    // printLog("temperature = %24.16g K\n", temperature[1][1][1]);

    // Calculate pressure.
    static gr_float *pressure;
    double pressure_units = grackle_code_units.density_units * pow(grackle_code_units.velocity_units, 2);
    if (once==0 && Dts!=NULL) {
        if (pressure!=NULL) FreeArray1D(pressure);
        pressure = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
    }
    if (Dts!=NULL) {
        if (calculate_pressure(&grackle_code_units, &grackle_chemistry_fields,
                            pressure) == 0) {
            printLog("call_grackle(): Error in calculate_pressure.\n");
            QUIT_PLUTO(1);
        }
    }
    // int count = 0;
    counter = 0;
    DOM_LOOP(k, j, i) {
        // if (count == 0)
        //     printLog("> step %d: lambda = %24.16e erg cm^3 s^-1\n", g_stepNumber, CONST_kB*d->Vgrac[TEMP][k][j][i]/((d->Vc[RHO][k][j][i]*UNIT_DENSITY/(d->Vgrac[MU][k][j][i]*CONST_mp))*fabs(cool_time_min) * grackle_code_units.time_units)/(g_gamma-1) );
        // count++;
        if (one_cell==0)
            id = (k-grid->lbeg[KDIR]) * grid->np_int[JDIR] * grid->np_int[IDIR] + (j-grid->lbeg[JDIR]) * grid->np_int[IDIR] + (i-grid->lbeg[IDIR]);
        if (one_cell==1) {
            if (counter>0) continue;
            i = cell_i;
            j = cell_j;
            k = cell_k;
        }
        // printLog("> step %d t=%e, dt=%e, before: prs/kB = %e, temp = %e, mu=%f\n", g_stepNumber, g_time, g_dt, (d->Vc[PRS][k][j][i]*UNIT_DENSITY*pow(UNIT_VELOCITY,2))/CONST_kB, d->Vgrac[TEMP][k][j][i], d->Vgrac[MU][k][j][i]);
        if (Dts!=NULL) d->Vc[PRS][k][j][i] = pressure[id];
        d->Vgrac[TEMP][k][j][i] = temperature[id];
        if (one_cell==1)
            d->Vgrac[MU][k][j][i] = (d->Vc[RHO][k][j][i]*UNIT_DENSITY)/(d->Vc[PRS][k][j][i]*UNIT_DENSITY*pow(UNIT_VELOCITY, 2))*(CONST_kB/CONST_mp)*d->Vgrac[TEMP][k][j][i];

        if (grackle_config_data->primordial_chemistry >= 1) {
            d->Vc[X_HI][k][j][i]    = grackle_chemistry_fields.HI_density[id]/(grackle_config_data->HydrogenFractionByMass * grackle_chemistry_fields.density[id]);
            d->Vc[X_HII][k][j][i]   = grackle_chemistry_fields.HII_density[id]/(grackle_chemistry_fields.density[id] * grackle_config_data->HydrogenFractionByMass);
            d->Vc[Y_HeI][k][j][i]   = grackle_chemistry_fields.HeI_density[id]/(grackle_chemistry_fields.density[id] * (1-grackle_config_data->HydrogenFractionByMass));
            d->Vc[Y_HeII][k][j][i]  = grackle_chemistry_fields.HeII_density[id]/(grackle_chemistry_fields.density[id] * (1-grackle_config_data->HydrogenFractionByMass));
            d->Vc[Y_HeIII][k][j][i] = grackle_chemistry_fields.HeIII_density[id]/(grackle_chemistry_fields.density[id] * (1-grackle_config_data->HydrogenFractionByMass));
            d->Vc[elec][k][j][i]    = grackle_chemistry_fields.e_density[id]*(CONST_me/CONST_mp);
        }
        if (grackle_config_data->primordial_chemistry >= 2) {
            d->Vc[X_HM][k][j][i]      =  grackle_chemistry_fields.HM_density[id]/(grackle_chemistry_fields.density[id] * grackle_config_data->HydrogenFractionByMass);
	        d->Vc[X_H2I][k][j][i]     =  grackle_chemistry_fields.H2I_density[id]/(grackle_chemistry_fields.density[id] * grackle_config_data->HydrogenFractionByMass);
	        d->Vc[X_H2II][k][j][i]    =  grackle_chemistry_fields.H2II_density[id]/(grackle_chemistry_fields.density[id] * grackle_config_data->HydrogenFractionByMass);
        }
        if (grackle_config_data->primordial_chemistry >= 3) {
            d->Vc[X_DI][k][j][i]     =  grackle_chemistry_fields.DI_density[id]/(grackle_chemistry_fields.density[id] * grackle_config_data->HydrogenFractionByMass);
	        d->Vc[X_DII][k][j][i]    =  grackle_chemistry_fields.DII_density[id]/(grackle_chemistry_fields.density[id] * grackle_config_data->HydrogenFractionByMass);
	        d->Vc[X_HDI][k][j][i]    =  grackle_chemistry_fields.HDI_density[id]/(grackle_chemistry_fields.density[id] * grackle_config_data->HydrogenFractionByMass);
        }
        normalize_ions_grackle(d, grackle_config_data, i, j, k);
        // solar metallicity
        if (grackle_config_data->metal_cooling == 1) 
            d->Vc[Z_MET][k][j][i] = grackle_chemistry_fields.metal_density[id]/(grackle_chemistry_fields.density[id]*grackle_config_data->SolarMetalFractionByMass);

        // printLog("> step %d t=%e, dt=%e,  after: prs/kB = %e, temp = %e, mu=%f\n", g_stepNumber, g_time, g_dt, (d->Vc[PRS][k][j][i]*UNIT_DENSITY*pow(UNIT_VELOCITY,2))/CONST_kB, d->Vgrac[TEMP][k][j][i], d->Vgrac[MU][k][j][i]);
        counter++;
    }
    if (one_cell==0)
        MeanMolecularWeight(d, grid);
    // printLog(stdout, "pressure = %24.16g dyne/cm^2\n", pressure[1][1][1]*pressure_units);

    /*
    // Calculate gamma.
    static gr_float *gamma;
    if (once==0)
        gamma = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
    if (calculate_gamma(&grackle_code_units, &grackle_chemistry_fields,
                        gamma) == 0) {
        printLog("Error in calculate_gamma.\n");
        QUIT_PLUTO(1);
    }

    // printLog("gamma = %24.16g\n", gamma[1][1][1]);
    
    // Calculate dust temperature.
    static gr_float *dust_temperature;
    if (once==0)
        dust_temperature = ARRAY_1D((one_cell==0)?(grid->np_int[KDIR] * grid->np_int[JDIR] * grid->np_int[IDIR]):1, gr_float);
    if (calculate_dust_temperature(&grackle_code_units, &grackle_chemistry_fields,
                      dust_temperature) == 0) {
        printLog("Error in calculate_dust_temperature.\n");
        QUIT_PLUTO(1);
    }

    // printLog("dust_temperature = %24.16g K\n", dust_temperature[1][1][1]);
    */
    if (one_cell==0) {
        once ++;
    }
    else {
        free_chemistry_data();
        once = 0;
    }
}
