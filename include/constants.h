/*
 * Public header file containing the constants
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <petsc.h>
#define freq 60.0                   /**< System frequency */
#define w_s (2.0 * PETSC_PI * freq) /**< Angular speed */
#define MAXLINE 10000               /**< Max. number of characters in a line */
#define ISOLATED_BUS 4              /**< Isolated bus */
#define REF_BUS 3                   /**< Reference bus (swing bus) */
#define PV_BUS 2                    /**< PV (voltage-controlled) bus */
#define PQ_BUS 1                    /**< PQ bus */
#define NGEN_AT_BUS_MAX                                                        \
  32                        /**< Maximum number of generators allowed at a bus \
                             */
#define NLOAD_AT_BUS_MAX 32 /**< Maximum number of loads allowed at a bus */
#define NPHASE 1            /**< Per-phase analysis */
#define MAX_KV_LEVELS 200   /**< Maximum KV levels */

/* Fuel types for generators */
#define GENFUEL_COAL 0      /* Coal */
#define GENFUEL_WIND 1      /* Wind */
#define GENFUEL_SOLAR 2     /* Solar */
#define GENFUEL_NG 3        /* Natural gas */
#define GENFUEL_NUCLEAR 4   /* Nuclear */
#define GENFUEL_HYDRO 5     /* hydro */
#define GENFUEL_UNDEFINED 6 /* All other sources (undefined) */

/* Default ramp rates (MW/min) */
/* From
 * https://www.researchgate.net/post/What_is_the_typical_MW_minute_ramping_capability_for_each_type_of_reserve
 */
#define GENRAMPRATE_COAL 3      /* Coal */
#define GENRAMPRATE_WIND 45     /* Wind */
#define GENRAMPRATE_SOLAR 200   /* Solar */
#define GENRAMPRATE_NG 35       /* Natural gas */
#define GENRAMPRATE_NUCLEAR 20  /* Nuclear */
#define GENRAMPRATE_HYDRO 150   /* Hydro */
#define GENRAMPRATE_UNDEFINED 3 /* All other sources (same as coal) */

#define EXAGO_IGNORE -1000000 /* Ignore value */

#define BAD_INPUT_DATA "Bad Input Data :"
#define INCORRECT_BOUNDS "Incorrect Bounds :"

/*
 * We often want a char[] with enough space to handle any solver or model name.
 * These sizes are greater than the greatest of the model or solver name,
 * respectively.
 */
static constexpr std::size_t max_model_name_len = 64;
static constexpr std::size_t max_solver_name_len = 64;

#endif
