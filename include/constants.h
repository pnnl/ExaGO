/*
 * Public header file containing the constants
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <petsc.h>
#define freq 60.0                   /**< System frequency */
#define w_s (2.0 * PETSC_PI * freq) /**< Angular speed */
#define MAXLINE 4095                /**< Max. number of characters in a line */
#define ISOLATED_BUS 4              /**< Isolated bus */
#define REF_BUS 3                   /**< Reference bus (swing bus) */
#define PV_BUS 2                    /**< PV (voltage-controlled) bus */
#define PQ_BUS 1                    /**< PQ bus */
#define NGEN_AT_BUS_MAX                                                        \
  32                        /**< Maximum number of generators allowed at a bus \
                             */
#define NLOAD_AT_BUS_MAX 32 /**< Maximum number of loads allowed at a bus */
#define NPHASE 1            /**< Per-phase analysis */

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

#endif
