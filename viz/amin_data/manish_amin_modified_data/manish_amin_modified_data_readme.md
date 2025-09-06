# Complete File-by-File Analysis: Manish Amin Modified Data

## Dataset Overview
**Total Files**: 711 files organized in comprehensive energy system analysis  
**Coverage**: Western Electricity Coordinating Council (WECC) region  
**Scenarios**: 4 case studies + 4 additional result sets  
**Time Resolution**: Hourly data (24 hours)  
**Geographic Scope**: 31 Balancing Authorities

## COMPLETE DIRECTORY TREE STRUCTURE (from viz folder)

```
viz/amin_data/manish_amin_modified_data/
├── WECC_BA_CAPACITY.csv
├── Total Opeation Cost.xlsx
├── manish_amin_modified_data_readme.md
├── Locations of Power Plants/
│   ├── Western_Power_plants_Locations-USA.csv
│   └── Western_Power_plants_Locations-CANADA-MEX_Update.csv
├── addtiional/
│   ├── results/
│   │   ├── data_center/
│   │   │   ├── dc_flexibility_hourly_by_ba.csv
│   │   │   ├── dc_flexibility_resource_detailed.csv
│   │   │   ├── dc_case1_load_after_flexibility.csv
│   │   │   └── dc_base_load.csv
│   │   ├── storage/
│   │   │   └── storage_operation.csv
│   │   ├── figures/
│   │   │   ├── battery_charge_discharge_hourly.png
│   │   │   ├── cost_breakdown_stacked.png
│   │   │   ├── generation_by_fuel_type_lines.png
│   │   │   ├── generation_by_fuel_type_stacked.png
│   │   │   ├── hourly_cost_comparison.png
│   │   │   └── load_vs_generation.png
│   │   ├── lmp_data/
│   │   │   └── lmp_by_ba_hour.csv
│   │   ├── transmission/
│   │   │   ├── ba_import_export_analysis.csv
│   │   │   ├── ba_import_export_summary.csv
│   │   │   ├── ba_import_export_violations_5000MW.csv
│   │   │   └── transmission_flows.csv
│   │   ├── summary_statistics.csv
│   │   ├── power_generation/
│   │   │   ├── AESO_generation_by_fuel.csv
│   │   │   ├── AVA_generation_by_fuel.csv
│   │   │   ├── AZPS_generation_by_fuel.csv
│   │   │   ├── BANC_generation_by_fuel.csv
│   │   │   ├── BCHA_generation_by_fuel.csv
│   │   │   ├── BPAT_generation_by_fuel.csv
│   │   │   ├── CENACE_generation_by_fuel.csv
│   │   │   ├── CHPD_generation_by_fuel.csv
│   │   │   ├── CISO_generation_by_fuel.csv
│   │   │   ├── DOPD_generation_by_fuel.csv
│   │   │   ├── EPE_generation_by_fuel.csv
│   │   │   ├── GCPD_generation_by_fuel.csv
│   │   │   ├── IID_generation_by_fuel.csv
│   │   │   ├── IPCO_generation_by_fuel.csv
│   │   │   ├── LDWP_generation_by_fuel.csv
│   │   │   ├── NEVP_generation_by_fuel.csv
│   │   │   ├── NWMT_generation_by_fuel.csv
│   │   │   ├── PACE_generation_by_fuel.csv
│   │   │   ├── PACW_generation_by_fuel.csv
│   │   │   ├── PGE_generation_by_fuel.csv
│   │   │   ├── PNM_generation_by_fuel.csv
│   │   │   ├── PSCO_generation_by_fuel.csv
│   │   │   ├── PSEI_generation_by_fuel.csv
│   │   │   ├── SCL_generation_by_fuel.csv
│   │   │   ├── SRP_generation_by_fuel.csv
│   │   │   ├── TEPC_generation_by_fuel.csv
│   │   │   ├── TIDC_generation_by_fuel.csv
│   │   │   ├── TPWR_generation_by_fuel.csv
│   │   │   ├── WACM_generation_by_fuel.csv
│   │   │   ├── WALC_generation_by_fuel.csv
│   │   │   ├── WAUW_generation_by_fuel.csv
│   │   │   ├── system_generation_by_fuel.csv
│   │   │   └── generation_schedule.csv
│   │   └── operation_costs/
│   │       ├── AESO_hourly_operation_costs.csv
│   │       ├── AVA_hourly_operation_costs.csv
│   │       ├── AZPS_hourly_operation_costs.csv
│   │       ├── BANC_hourly_operation_costs.csv
│   │       ├── BCHA_hourly_operation_costs.csv
│   │       ├── BPAT_hourly_operation_costs.csv
│   │       ├── CENACE_hourly_operation_costs.csv
│   │       ├── CHPD_hourly_operation_costs.csv
│   │       ├── CISO_hourly_operation_costs.csv
│   │       ├── DOPD_hourly_operation_costs.csv
│   │       ├── EPE_hourly_operation_costs.csv
│   │       ├── GCPD_hourly_operation_costs.csv
│   │       ├── IID_hourly_operation_costs.csv
│   │       ├── IPCO_hourly_operation_costs.csv
│   │       ├── LDWP_hourly_operation_costs.csv
│   │       ├── NEVP_hourly_operation_costs.csv
│   │       ├── NWMT_hourly_operation_costs.csv
│   │       ├── PACE_hourly_operation_costs.csv
│   │       ├── PACW_hourly_operation_costs.csv
│   │       ├── PGE_hourly_operation_costs.csv
│   │       ├── PNM_hourly_operation_costs.csv
│   │       ├── PSCO_hourly_operation_costs.csv
│   │       ├── PSEI_hourly_operation_costs.csv
│   │       ├── SCL_hourly_operation_costs.csv
│   │       ├── SRP_hourly_operation_costs.csv
│   │       ├── TEPC_hourly_operation_costs.csv
│   │       ├── TIDC_hourly_operation_costs.csv
│   │       ├── TPWR_hourly_operation_costs.csv
│   │       ├── WACM_hourly_operation_costs.csv
│   │       ├── WALC_hourly_operation_costs.csv
│   │       ├── WAUW_hourly_operation_costs.csv
│   │       ├── ba_cost_summary.csv
│   │       ├── generator_operation_costs.csv
│   │       └── hourly_system_operation_costs.csv
│   ├── results_2/
│   │   └── [IDENTICAL STRUCTURE TO results/]
│   ├── results_3/
│   │   └── [IDENTICAL STRUCTURE TO results/]
│   └── results_4/
│       └── [IDENTICAL STRUCTURE TO results/]
├── Case study_1/
│   └── CASE STUDY 1/
│       ├── Balancing Authority Demand (MW).csv
│       ├── Data Center Demand (MW).csv
│       ├── Balancing Authority Hourly Operation Costs/
│       │   ├── [31 BA files: {BA}_hourly_operation_costs.csv]
│       │   ├── ba_cost_summary.csv
│       │   ├── generator_operation_costs.csv
│       │   └── hourly_system_operation_costs.csv
│       ├── Balancing Authority Power Generation/
│       │   ├── [31 BA files: {BA}_generation_by_fuel.csv]
│       │   ├── system_generation_by_fuel.csv
│       │   └── generation_schedule.csv
│       ├── Zonal Price/
│       │   └── lmp_by_ba_hour.csv
│       └── Additional files/
│           ├── Balancing Authority Hourly Operation Costs - Copy/
│           │   └── [34 cost files - copies]
│           ├── figures/
│           │   └── [6 PNG visualization files]
│           ├── storage/
│           │   └── storage_operation.csv
│           ├── summary_statistics.csv
│           └── transmission/
│               └── [4 transmission analysis files]
├── Case study_2/
│   ├── Balancing Authority Demand (MW).csv
│   ├── Data Center Demand (MW).csv
│   ├── Balancing Authority Hourly Operation Costs/
│   │   └── [34 cost files]
│   ├── Balancing Authority Power Generation/
│   │   └── [33 generation files]
│   ├── Data Center Flexibility/
│   │   └── Data Center Energy Flexibility.csv
│   ├── Zonal Price/
│   │   └── lmp_by_ba_hour.csv
│   └── Additional/
│       ├── dc_base_load.csv
│       ├── dc_case1_load_after_flexibility.csv
│       ├── dc_flexibility_resource_detailed.csv
│       ├── figures/
│       │   └── [6 PNG files]
│       ├── storage/
│       │   └── storage_operation.csv
│       ├── summary_statistics.csv
│       └── transmission/
│           └── [4 transmission files]
├── Case study_3/
│   └── [IDENTICAL STRUCTURE TO Case study_2/]
└── Case study_4/
    └── [IDENTICAL STRUCTURE TO Case study_2/]
```

### File Path Templates for Programming

**Root Level Files:**
```
viz/amin_data/manish_amin_modified_data/WECC_BA_CAPACITY.csv
viz/amin_data/manish_amin_modified_data/Total Opeation Cost.xlsx
viz/amin_data/manish_amin_modified_data/Locations of Power Plants/Western_Power_plants_Locations-USA.csv
viz/amin_data/manish_amin_modified_data/Locations of Power Plants/Western_Power_plants_Locations-CANADA-MEX_Update.csv
```

**Additional Results (4 identical sets):**
```
viz/amin_data/manish_amin_modified_data/addtiional/results/[subfolder]/[filename]
viz/amin_data/manish_amin_modified_data/addtiional/results_2/[subfolder]/[filename]
viz/amin_data/manish_amin_modified_data/addtiional/results_3/[subfolder]/[filename]
viz/amin_data/manish_amin_modified_data/addtiional/results_4/[subfolder]/[filename]
```

**Case Studies:**
```
viz/amin_data/manish_amin_modified_data/Case study_1/CASE STUDY 1/[subfolder]/[filename]
viz/amin_data/manish_amin_modified_data/Case study_2/[subfolder]/[filename]
viz/amin_data/manish_amin_modified_data/Case study_3/[subfolder]/[filename]
viz/amin_data/manish_amin_modified_data/Case study_4/[subfolder]/[filename]
```

**Balancing Authority File Patterns:**
```
# Individual BA generation files:
viz/amin_data/manish_amin_modified_data/[scenario]/[subfolder]/power_generation/{BA}_generation_by_fuel.csv

# Individual BA cost files:
viz/amin_data/manish_amin_modified_data/[scenario]/[subfolder]/operation_costs/{BA}_hourly_operation_costs.csv

# Where {BA} = AESO, AVA, AZPS, BANC, BCHA, BPAT, CENACE, CHPD, CISO, DOPD, EPE, GCPD, IID, IPCO, LDWP, NEVP, NWMT, PACE, PACW, PGE, PNM, PSCO, PSEI, SCL, SRP, TEPC, TIDC, TPWR, WACM, WALC, WAUW
```

---

## ROOT LEVEL FILES (3 files)

### 1. `WECC_BA_CAPACITY.csv`
- **Purpose**: Master capacity reference for all Balancing Authorities
- **Structure**: 32 rows × 10 columns
- **Columns**: BA, Hydro, Nuclear, Coal, Natural Gas, Geothermal, Biomass, Wind, PV, Battery_Storage
- **Data Range**: 0 to 31,686 MW capacity values
- **Usage**: Foundation data for all capacity constraints

### 2. `Total Opeation Cost.xlsx`
- **Purpose**: Excel summary of operation costs across scenarios
- **Format**: Microsoft Excel workbook
- **Content**: Aggregated cost comparisons between case studies
- **Usage**: High-level financial analysis

### 3. `manish_amin_modified_data_readme.md`
- **Purpose**: This comprehensive documentation file
- **Content**: Analysis of all 711 files in dataset

---

## LOCATIONS OF POWER PLANTS FOLDER (2 files)

### 4. `Western_Power_plants_Locations-USA.csv`
- **Purpose**: Geographic and technical data for US power plants
- **Structure**: 3,349 plants × 18 columns
- **Key Columns**: Plant Code, Name, Lat/Lon, State, County, BA, Primary Type, Total Capacity, 9 fuel-specific capacities
- **Coverage**: All western US states in WECC region
- **Usage**: Plant mapping, capacity visualization, geographic analysis

### 5. `Western_Power_plants_Locations-CANADA-MEX_Update.csv`
- **Purpose**: International power plants in WECC region
- **Structure**: 277 plants × 10 columns
- **Key Addition**: Country field (CAN/MEX)
- **Coverage**: Cross-border plants connected to WECC
- **Usage**: International power flow analysis

---

## ADDITIONAL RESULTS FOLDER (404 files total)

Contains 4 identical result sets: `results/`, `results_2/`, `results_3/`, `results_4/`  
Each set has 101 files with identical structure but different scenario parameters.

### DATA CENTER SUBFOLDER (4 files per result set = 16 files)

#### Files 6-9: `dc_flexibility_hourly_by_ba.csv` (×4 copies)
- **Purpose**: Hourly data center flexibility by BA
- **Structure**: 24 hours × 13 columns (Hour + 12 BA flexibility values)
- **Data Range**: -880.57 to +440.29 MW
- **Interpretation**: Negative = demand increase, Positive = demand reduction
- **Usage**: Data center demand response analysis

#### Files 10-13: `dc_flexibility_resource_detailed.csv` (×4 copies)
- **Purpose**: Detailed flexibility resource breakdown
- **Structure**: 288 rows × 6 columns
- **Coverage**: 12 BAs × 24 hours
- **Key Columns**: Hour, BA, Flexibility_Resource_MW, Base_Total_Load_MW, Net_Load_MW, Flexibility_Type
- **Usage**: Granular flexibility tracking

#### Files 14-17: `dc_case1_load_after_flexibility.csv` (×4 copies)
- **Purpose**: Final data center loads post-flexibility
- **Structure**: 288 rows × 7 columns
- **Key Metrics**: Final_Server_Load_MW, Final_Cooling_Load_MW, Server_Capacity_Utilization_%
- **Usage**: Post-flexibility impact analysis

#### Files 18-21: `dc_base_load.csv` (×4 copies)
- **Purpose**: Baseline data center loads pre-flexibility
- **Structure**: 288 rows × 8 columns
- **Key Metrics**: DC_Grid_Capacity_MW, Server_Max_Capacity_MW, Workload_Pattern, Base loads
- **Usage**: Baseline comparison for flexibility studies

### STORAGE SUBFOLDER (4 files)

#### Files 22-25: `storage_operation.csv` (×4 copies)
- **Purpose**: Battery storage operations by BA and hour
- **Structure**: 216 rows × 6 columns (9 BAs with storage × 24 hours)
- **Key Columns**: Hour, BA, Charge_MW, Discharge_MW, SOC_MWh, Capacity_MWh
- **Data Range**: Charge/Discharge 0-3,371 MW, Capacity 352-29,114 MWh
- **Usage**: Storage optimization analysis

### FIGURES SUBFOLDER (24 PNG files)

#### Files 26-49: Visualization files (6 per result set)
- `battery_charge_discharge_hourly.png`: Storage operation patterns
- `cost_breakdown_stacked.png`: Cost component analysis
- `generation_by_fuel_type_lines.png`: Generation trends by fuel
- `generation_by_fuel_type_stacked.png`: Generation mix composition
- `hourly_cost_comparison.png`: Hourly cost variations
- `load_vs_generation.png`: System load vs generation balance

### LMP DATA SUBFOLDER (4 files)

#### Files 50-53: `lmp_by_ba_hour.csv` (×4 copies)
- **Purpose**: Locational Marginal Prices by BA and hour
- **Structure**: 24 hours × 32 columns (Hour + 31 BA LMP values)
- **Price Range**: $24.82 - $317.52/MWh
- **Peak Prices**: Hours 10-16 (high demand periods)
- **Usage**: Economic dispatch verification, price analysis

### TRANSMISSION SUBFOLDER (16 files)

#### Files 54-57: `ba_import_export_analysis.csv` (×4 copies)
- **Purpose**: Detailed import/export analysis by BA
- **Structure**: 744 rows × 13 columns (31 BAs × 24 hours)
- **Key Columns**: Load_MW, Local_Supply_MW, Import_Export_MW, Status, Exceeds_5000MW_Limit
- **Usage**: Transmission constraint analysis

#### Files 58-61: `ba_import_export_summary.csv` (×4 copies)
- **Purpose**: Summary of import/export flows
- **Structure**: Variable rows × 6 columns
- **Key Columns**: Hour, From_BA, To_BA, Flow_MW, Capacity_MW, Utilization_%
- **Usage**: Transmission utilization analysis

#### Files 62-65: `ba_import_export_violations_5000MW.csv` (×4 copies)
- **Purpose**: Tracking 5000MW import limit violations
- **Structure**: Variable rows × 8 columns (only violating BAs)
- **Key Focus**: CISO frequently exceeds 5000MW import limit
- **Usage**: Transmission constraint violation analysis

#### Files 66-69: `transmission_flows.csv` (×4 copies)
- **Purpose**: Detailed transmission flows between all BAs
- **Structure**: ~44,400 rows × 6 columns (~1,850 lines × 24 hours)
- **Coverage**: All transmission paths between BAs
- **Usage**: Detailed power flow analysis

### SUMMARY STATISTICS SUBFOLDER (4 files)

#### Files 70-73: `summary_statistics.csv` (×4 copies)
- **Purpose**: High-level system performance metrics
- **Structure**: 1 row × 36 columns (daily summary)
- **Key Metrics**: 
  - Total System Cost: ~$52.8M/day
  - Load Shedding: 0 MWh (system reliable)
  - Renewable Curtailment: 0 MWh (full utilization)
- **Usage**: System performance benchmarking

### POWER GENERATION SUBFOLDER (132 files)

#### Files 74-205: Individual BA generation files (31 BAs × 4 result sets)
**File Pattern**: `[BA]_generation_by_fuel.csv`
**BA List**: AESO, AVA, AZPS, BANC, BCHA, BPAT, CENACE, CHPD, CISO, DOPD, EPE, GCPD, IID, IPCO, LDWP, NEVP, NWMT, PACE, PACW, PGE, PNM, PSCO, PSEI, SCL, SRP, TEPC, TIDC, TPWR, WACM, WALC, WAUW

- **Purpose**: Hourly generation by fuel type for each BA
- **Structure**: 24 hours × 11 columns per file
- **Columns**: Hour, NG_MW, GEO_MW, BIO_MW, NUCLEAR_MW, COAL_MW, WIND_MW, PV_MW, HYDRO_MW, BATTERY_MW, IMPORT/EXPORT_MW
- **Usage**: Individual BA fuel mix analysis

#### Files 206-209: `system_generation_by_fuel.csv` (×4 copies)
- **Purpose**: System-wide generation totals by fuel type
- **Structure**: 24 hours × 11 columns (same as individual files)
- **Usage**: System-level generation analysis

#### Files 210-213: `generation_schedule.csv` (×4 copies)
- **Purpose**: Individual generator dispatch schedule
- **Structure**: ~17,856 rows × 7 columns
- **Key Columns**: Hour, GeneratorKey, BA, Fuel_Type, Power_MW, Status, Startup
- **Usage**: Unit commitment analysis, generator-level dispatch

### OPERATION COSTS SUBFOLDER (136 files)

#### Files 214-349: Individual BA cost files (31 BAs × 4 result sets)
**File Pattern**: `[BA]_hourly_operation_costs.csv`

- **Purpose**: Detailed hourly cost breakdown by BA
- **Structure**: 24 hours × 7 columns per file
- **Key Columns**: BA_Startup_Costs_$, BA_Minimum_Fuel_Costs_$, BA_Variable_Costs_$, BA_Load_Shedding_Costs_$, BA_Import_Export_Costs_$, BA_Total_Hourly_Cost_$
- **Cost Range**: -$3M to +$820K per BA per hour
- **Usage**: Regional economic analysis

#### Files 350-353: `ba_cost_summary.csv` (×4 copies)
- **Purpose**: Daily cost summary by BA
- **Structure**: 31 BAs × 7 columns
- **Usage**: BA-level cost comparison

#### Files 354-357: `generator_operation_costs.csv` (×4 copies)
- **Purpose**: Individual generator cost details
- **Structure**: ~17,856 rows × 12 columns
- **Key Columns**: GeneratorKey, Power_MW, Startup costs, Fuel costs, Variable costs
- **Usage**: Unit-level economic analysis

#### Files 358-361: `hourly_system_operation_costs.csv` (×4 copies)
- **Purpose**: System-wide hourly cost totals
- **Structure**: 24 hours × 7 columns
- **Usage**: System-level economic analysis

---

## CASE STUDY 1 FOLDER (80 files)

### Main Files (2 files)

#### File 362: `Balancing Authority Demand (MW).csv`
- **Purpose**: Hourly demand by BA for baseline case
- **Structure**: 24 hours × 36 columns (Year, Month, Day, Period + 32 BA demand columns)
- **Data Range**: 100 - 35,000+ MW per BA per hour
- **Date**: July 16, 2026
- **Usage**: Base case demand analysis

#### File 363: `Data Center Demand (MW).csv`
- **Purpose**: Data center load patterns without flexibility
- **Structure**: 288 rows × 8 columns (12 BAs × 24 hours)
- **Usage**: Baseline data center analysis

### Balancing Authority Hourly Operation Costs Subfolder (34 files)

#### Files 364-397: Individual BA hourly cost files (31 BAs + 3 summary files)
Same structure as additional results operation costs

### Balancing Authority Power Generation Subfolder (33 files)

#### Files 398-430: Individual BA generation files (31 BAs + 2 system files)
Same structure as additional results power generation

### Zonal Price Subfolder (1 file)

#### File 431: `lmp_by_ba_hour.csv`
Same structure as additional results LMP files

### Additional Files Subfolder (46 files)

#### Files 432-477: Duplicate analysis files
Contains copies of operation costs, transmission data, storage data, summary statistics, and visualization figures

---

## CASE STUDY 2 FOLDER (73 files)

### Main Files (2 files)

#### File 478: `Balancing Authority Demand (MW).csv`
Case 2 specific demand patterns

#### File 479: `Data Center Demand (MW).csv`
Case 2 data center loads

### Data Center Flexibility Subfolder (1 file)

#### File 480: `Data Center Energy Flexibility.csv`
- **Purpose**: Data center flexibility deployment in Case 2
- **Structure**: Multiple flexibility metrics and scenarios
- **Usage**: Flexibility resource analysis

### Other Subfolders (70 files)
- Balancing Authority Hourly Operation Costs/ (34 files)
- Balancing Authority Power Generation/ (33 files)
- Zonal Price/ (1 file)
- Additional/ (15 files including storage, transmission, figures)

---

## CASE STUDY 3 FOLDER (73 files)

#### Files 481-553: Same structure as Case Study 2
Different scenario parameters but identical file organization

---

## CASE STUDY 4 FOLDER (73 files)

#### Files 554-626: Same structure as Case Studies 2 & 3
Different scenario parameters but identical file organization

---

## COMPREHENSIVE FILE SUMMARY

### File Type Distribution
- **CSV Files**: 663 (data tables)
- **PNG Files**: 48 (visualizations)
- **Excel Files**: 1 (summary)
- **Markdown Files**: 1 (documentation)
- **TOTAL**: 711 files

### Data Category Breakdown
1. **Capacity Data**: 1 file (WECC capacity reference)
2. **Demand Data**: 8 files (load patterns across scenarios)
3. **Generation Data**: 165 files (power output by fuel and BA)
4. **Cost Data**: 169 files (economic analysis at multiple levels)
5. **Price Data**: 8 files (LMP analysis across scenarios)
6. **Transmission Data**: 32 files (power flows and constraints)
7. **Storage Data**: 8 files (battery operation patterns)
8. **Data Center Data**: 32 files (demand response analysis)
9. **Summary Data**: 8 files (high-level system metrics)
10. **Visualization Data**: 48 files (chart images)
11. **Location Data**: 2 files (plant geographic information)

### Key Balancing Authorities (31 total)
AESO (Alberta), AVA (Avista), AZPS (Arizona Public Service), BANC (Balancing Authority of Northern California), BCHA (BC Hydro), BPAT (Bonneville Power Administration), CENACE (Mexico), CHPD (Chelan County PUD), CISO (California ISO), DOPD (Douglas County PUD), EPE (El Paso Electric), GCPD (Grant County PUD), IID (Imperial Irrigation District), IPCO (Idaho Power), LDWP (Los Angeles Department of Water and Power), NEVP (Nevada Power), NWMT (NorthWestern Energy Montana), PACE (PacifiCorp East), PACW (PacifiCorp West), PGE (Portland General Electric), PNM (Public Service Company of New Mexico), PSCO (Public Service Company of Colorado), PSEI (Puget Sound Energy), SCL (Seattle City Light), SRP (Salt River Project), TEPC (Tucson Electric Power), TIDC (Turlock Irrigation District), TPWR (Tacoma Power), WACM (Western Area Colorado Missouri), WALC (Western Area Lower Colorado), WAUW (Western Area Upper Great Plains West)

---

## DASHBOARD DEVELOPMENT GUIDELINES

### Recommended File Paths
```javascript
// Primary data sources
const capacityData = "WECC_BA_CAPACITY.csv"
const plantLocations = "Locations of Power Plants/"
const primaryScenario = "Case study_2/" // Main analysis
const detailedResults = "addtiional/results_2/" // Supplementary analysis
```

### Data Integration Points
- **BA codes**: Link all datasets (31 unique BAs)
- **Time stamps**: Hourly resolution (1-24)
- **Geographic data**: Lat/lon coordinates for mapping
- **Economic data**: Costs in USD, prices in $/MWh
- **Energy data**: Power in MW, Energy in MWh

### Visualization Opportunities
1. **Geographic Maps**: Plant locations with capacity bubbles
2. **Time Series**: Generation, demand, prices by hour
3. **Economic Analysis**: Cost breakdowns, LMP heatmaps
4. **System Operations**: Storage cycles, transmission flows
5. **Scenario Comparison**: Side-by-side case study analysis

## QUICK REFERENCE: COMMON FILE PATHS FOR GRAPH DEVELOPMENT

### Essential Data Files for Dashboard/Graphs

**System Capacity Data:**
```javascript
const capacityFile = "viz/amin_data/manish_amin_modified_data/WECC_BA_CAPACITY.csv";
```

**Plant Location Data:**
```javascript
const usaPlantsFile = "viz/amin_data/manish_amin_modified_data/Locations of Power Plants/Western_Power_plants_Locations-USA.csv";
const intlPlantsFile = "viz/amin_data/manish_amin_modified_data/Locations of Power Plants/Western_Power_plants_Locations-CANADA-MEX_Update.csv";
```

**Primary Case Study Data (Case Study 2 - Recommended):**
```javascript
// Main demand and generation data
const demandFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Balancing Authority Demand (MW).csv";
const dcDemandFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Data Center Demand (MW).csv";
const dcFlexibilityFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Data Center Flexibility/Data Center Energy Flexibility.csv";

// Pricing data
const lmpFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Zonal Price/lmp_by_ba_hour.csv";

// Individual BA generation (replace {BA} with actual BA code)
const baGenerationFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Balancing Authority Power Generation/{BA}_generation_by_fuel.csv";

// Individual BA costs (replace {BA} with actual BA code)
const baCostFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Balancing Authority Hourly Operation Costs/{BA}_hourly_operation_costs.csv";

// System-wide data
const systemGenerationFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Balancing Authority Power Generation/system_generation_by_fuel.csv";
const costSummaryFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Balancing Authority Hourly Operation Costs/ba_cost_summary.csv";

// Storage and transmission
const storageFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Additional/storage/storage_operation.csv";
const transmissionFlowsFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Additional/transmission/transmission_flows.csv";

// Summary statistics
const summaryStatsFile = "viz/amin_data/manish_amin_modified_data/Case study_2/Additional/summary_statistics.csv";
```

**Alternative Results Data (More Detailed Analysis):**
```javascript
// Use results_2 for most detailed analysis (matches Case study_2)
const detailedResultsBase = "viz/amin_data/manish_amin_modified_data/addtiional/results_2/";

// Data center analysis
const dcFlexHourlyFile = detailedResultsBase + "data_center/dc_flexibility_hourly_by_ba.csv";
const dcFlexDetailedFile = detailedResultsBase + "data_center/dc_flexibility_resource_detailed.csv";
const dcBaseLoadFile = detailedResultsBase + "data_center/dc_base_load.csv";
const dcAfterFlexFile = detailedResultsBase + "data_center/dc_case1_load_after_flexibility.csv";

// Transmission analysis
const importExportAnalysisFile = detailedResultsBase + "transmission/ba_import_export_analysis.csv";
const transmissionViolationsFile = detailedResultsBase + "transmission/ba_import_export_violations_5000MW.csv";

// Generator-level data
const generatorScheduleFile = detailedResultsBase + "power_generation/generation_schedule.csv";
const generatorCostsFile = detailedResultsBase + "operation_costs/generator_operation_costs.csv";

// System-level hourly costs
const hourlySystemCostsFile = detailedResultsBase + "operation_costs/hourly_system_operation_costs.csv";
```

**Pre-generated Visualizations:**
```javascript
// Chart images (PNG files) - available in multiple locations
const chartsBase = "viz/amin_data/manish_amin_modified_data/Case study_2/Additional/figures/";
const charts = {
    batteryOperation: chartsBase + "battery_charge_discharge_hourly.png",
    costBreakdown: chartsBase + "cost_breakdown_stacked.png",
    generationLines: chartsBase + "generation_by_fuel_type_lines.png",
    generationStacked: chartsBase + "generation_by_fuel_type_stacked.png",
    hourlyCosts: chartsBase + "hourly_cost_comparison.png",
    loadVsGeneration: chartsBase + "load_vs_generation.png"
};
```

### Dynamic File Path Construction

**For Multiple BAs (Loop through all Balancing Authorities):**
```javascript
const balancingAuthorities = [
    'AESO', 'AVA', 'AZPS', 'BANC', 'BCHA', 'BPAT', 'CENACE', 'CHPD', 
    'CISO', 'DOPD', 'EPE', 'GCPD', 'IID', 'IPCO', 'LDWP', 'NEVP', 
    'NWMT', 'PACE', 'PACW', 'PGE', 'PNM', 'PSCO', 'PSEI', 'SCL', 
    'SRP', 'TEPC', 'TIDC', 'TPWR', 'WACM', 'WALC', 'WAUW'
];

// Generate file paths for all BAs
const generationFiles = balancingAuthorities.map(ba => 
    `viz/amin_data/manish_amin_modified_data/Case study_2/Balancing Authority Power Generation/${ba}_generation_by_fuel.csv`
);

const costFiles = balancingAuthorities.map(ba => 
    `viz/amin_data/manish_amin_modified_data/Case study_2/Balancing Authority Hourly Operation Costs/${ba}_hourly_operation_costs.csv`
);
```

**For Multiple Scenarios (Compare across case studies):**
```javascript
const scenarios = ['Case study_1/CASE STUDY 1', 'Case study_2', 'Case study_3', 'Case study_4'];

// Generate demand file paths for all scenarios
const demandFiles = scenarios.map(scenario => 
    `viz/amin_data/manish_amin_modified_data/${scenario}/Balancing Authority Demand (MW).csv`
);

// Generate LMP file paths for all scenarios
const lmpFiles = scenarios.map(scenario => 
    `viz/amin_data/manish_amin_modified_data/${scenario}/Zonal Price/lmp_by_ba_hour.csv`
);
```

**For Multiple Result Sets (Detailed Analysis):**
```javascript
const resultSets = ['results', 'results_2', 'results_3', 'results_4'];

// Generate summary statistics paths for all result sets
const summaryFiles = resultSets.map(resultSet => 
    `viz/amin_data/manish_amin_modified_data/addtiional/${resultSet}/summary_statistics.csv`
);
```

### File Structure Validation Helper

```javascript
// Helper function to construct and validate file paths
function getDataFilePath(category, scenario = 'Case study_2', ba = null, resultSet = 'results_2') {
    const basePath = 'viz/amin_data/manish_amin_modified_data/';
    
    switch(category) {
        case 'capacity':
            return basePath + 'WECC_BA_CAPACITY.csv';
        
        case 'plants_usa':
            return basePath + 'Locations of Power Plants/Western_Power_plants_Locations-USA.csv';
        
        case 'demand':
            return basePath + `${scenario}/Balancing Authority Demand (MW).csv`;
        
        case 'ba_generation':
            if (!ba) throw new Error('BA code required for ba_generation');
            return basePath + `${scenario}/Balancing Authority Power Generation/${ba}_generation_by_fuel.csv`;
        
        case 'ba_costs':
            if (!ba) throw new Error('BA code required for ba_costs');
            return basePath + `${scenario}/Balancing Authority Hourly Operation Costs/${ba}_hourly_operation_costs.csv`;
        
        case 'lmp':
            return basePath + `${scenario}/Zonal Price/lmp_by_ba_hour.csv`;
        
        case 'storage':
            return basePath + `addtiional/${resultSet}/storage/storage_operation.csv`;
        
        case 'transmission':
            return basePath + `addtiional/${resultSet}/transmission/transmission_flows.csv`;
        
        case 'summary':
            return basePath + `addtiional/${resultSet}/summary_statistics.csv`;
        
        default:
            throw new Error(`Unknown category: ${category}`);
    }
}

// Usage examples:
// const capacityPath = getDataFilePath('capacity');
// const cisoGenPath = getDataFilePath('ba_generation', 'Case study_2', 'CISO');
// const lmpPath = getDataFilePath('lmp', 'Case study_3');
```

---

**COMPLETE DATASET DOCUMENTED: All 711 files analyzed and catalogued**
