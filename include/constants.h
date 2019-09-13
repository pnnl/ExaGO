/**
 * @file constants.h
 * @brief Public header file defining constants
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <petsc.h>
#define freq 60.0 /**< System frequency */
#define w_s (2.0*PETSC_PI*freq) /**< Angular speed */
#define epsilon 1E-8 /**< epsilon value for finite differences */
#define MAXLINE 1000 /**< Max. number of characters in a line */
#define ISOLATED_BUS 4 /**< Isolated bus */
#define REF_BUS 3 /**< Reference bus (swing bus) */
#define PV_BUS 2 /**< PV (voltage-controlled) bus */
#define PQ_BUS 1 /**< PQ bus */
#define NGEN_AT_BUS_MAX 32 /**< Maximum number of generators allowed at a bus */
#define NLOAD_AT_BUS_MAX 32 /**< Maximum number of loads allowed at a bus */
#define NPHASE 1 /**< Per-phase analysis */

#endif
