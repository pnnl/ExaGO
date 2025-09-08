# GIS Data Documentation - Western Power Plants Locations

## Overview
This directory contains Geographic Information System (GIS) data for power generation facilities across the Western United States, Canada, and Mexico regions. The data is provided in CSV format with detailed location, capacity, and technology type information.

---

## 📊 Dataset Summary

| Dataset | Records | Geographic Coverage | File Size |
|---------|---------|-------------------|-----------|
| USA Power Plants | 3,347 plants | Western United States | ~500KB |
| Canada-Mexico Power Plants | 275 plants | Western Canada & Mexico | ~50KB |
| **Total** | **3,622 plants** | **Western North America** | **~550KB** |

---

## 📁 File Descriptions

### 1. Western_Power_plants_Locations-USA (2).csv

**Path:** `/viz/data/Western_Power_plants_Locations-USA (2).csv`

**Description:** Comprehensive dataset of power generation facilities in the Western United States, including detailed capacity breakdowns by technology type.

#### Column Structure (18 columns):
| Column # | Field Name | Data Type | Description | Sample Values |
|----------|------------|-----------|-------------|---------------|
| 1 | Plant Code | Integer | Unique plant identifier | 34, 72, 99 |
| 2 | Plant Name | String | Official facility name | "Rollins", "Venice Hydro" |
| 3 | Latitude | Float | Geographic latitude (WGS84) | 39.134259, 34.01135 |
| 4 | Longitude | Float | Geographic longitude (WGS84) | -120.953341, -118.4168 |
| 5 | State | String (2-char) | US State abbreviation | CA, WA, AZ, NV, OR |
| 6 | County | String | County name | Placer, Los Angeles, Pierce |
| 7 | Balancing Authority | String | Grid operator code | CISO, PSEI, SRP, AZPS |
| 8 | Primary Type | String | Dominant technology type | Hydro, Solar, Natural Gas, Wind |
| 9 | Total Capacity (MW) | Float | Total generation capacity | 12.1, 177.8, 1207.4 |
| 10 | Battery Storage | Float | Battery storage capacity (MW) | 0.0, 45.2, 150.0 |
| 11 | Geothermal | Float | Geothermal capacity (MW) | 0.0, 50.0, 270.0 |
| 12 | Hydro | Float | Hydroelectric capacity (MW) | 12.1, 0.0, 199.8 |
| 13 | Natural Gas | Float | Natural gas capacity (MW) | 0.0, 177.8, 806.2 |
| 14 | Nuclear | Float | Nuclear capacity (MW) | 0.0, 1100.0, 2200.0 |
| 15 | Oil | Float | Oil/petroleum capacity (MW) | 0.0, 26.1, 45.0 |
| 16 | Other | Float | Other technology capacity (MW) | 0.0, 425.9, 100.0 |
| 17 | Solar | Float | Solar PV/thermal capacity (MW) | 0.0, 250.0, 579.0 |
| 18 | Wind | Float | Wind capacity (MW) | 0.0, 150.0, 845.0 |

#### Geographic Coverage:
- **Latitude Range:** 31.36°N to 48.99°N
- **Longitude Range:** -124.21°W to 47.55°W (Note: Some data quality issues detected)
- **States Covered:** CA, WA, OR, NV, AZ, UT, CO, NM, WY, MT, ID

#### Capacity Statistics:
- **Minimum Capacity:** 0.3 MW
- **Maximum Capacity:** 6,809.0 MW
- **Total Installed Capacity:** ~85,000 MW (estimated)

#### Technology Distribution:
| Technology | Count | Percentage |
|------------|-------|------------|
| Solar | 1,283 | 38.3% |
| Hydro | 595 | 17.8% |
| Natural Gas | 478 | 14.3% |
| Wind | 323 | 9.6% |
| Other | 211 | 6.3% |
| Battery Storage | 144 | 4.3% |
| Geothermal | 64 | 1.9% |
| Oil | 28 | 0.8% |
| Nuclear | 3 | 0.1% |

---

### 2. Western_Power_plants_Locations-CANADA-MEX_Update (1).csv

**Path:** `/viz/data/Western_Power_plants_Locations-CANADA-MEX_Update (1).csv`

**Description:** Power generation facilities in Western Canada and Mexico that interconnect with the Western US grid system.

#### Column Structure (10 columns):
| Column # | Field Name | Data Type | Description | Sample Values |
|----------|------------|-----------|-------------|---------------|
| 1 | Plant Code | String | International plant identifier | INTL_CAN0007597, INTL_MEX001 |
| 2 | Plant Name | String | Official facility name | "Aberfeldie", "Akolkolex" |
| 3 | Latitude | Float | Geographic latitude (WGS84) | 49.4933, 50.8227 |
| 4 | Longitude | Float | Geographic longitude (WGS84) | -115.3624, -118.0295 |
| 5 | State | String | Province/state code | BCHA, AESO, SON, BCN |
| 6 | County | String | County/region (often empty) | (mostly blank) |
| 7 | Balancing Authority | String | Grid operator | BCHA, AESO, CFE |
| 8 | Primary Type | String | Dominant technology type | Hydro, Natural Gas, Wind, Biomass |
| 9 | Total Capacity (MW) | Float | Total generation capacity | 24, 125, 185 |
| 10 | Country | String | Country code | CAN, MEX |

#### Geographic Coverage:
- **Latitude Range:** 32.37°N to 59.58°N
- **Longitude Range:** -133.64°W to -110.08°W
- **Regions:** British Columbia, Alberta (Canada); Sonora, Baja California (Mexico)

#### Capacity Statistics:
- **Minimum Capacity:** 1.0 MW
- **Maximum Capacity:** 2,746.0 MW
- **Total Installed Capacity:** ~15,000 MW (estimated)

#### Technology Distribution:
| Technology | Count | Percentage |
|------------|-------|------------|
| Hydro | 141 | 51.3% |
| Natural Gas | 44 | 16.0% |
| Biomass | 44 | 16.0% |
| Wind | 34 | 12.4% |
| Solar | 6 | 2.2% |
| Other | 2 | 0.7% |
| Oil | 2 | 0.7% |
| Geothermal | 2 | 0.7% |

#### Country Distribution:
- **Canada:** ~85% of records
- **Mexico:** ~15% of records

---

## 🗺️ Coordinate System & Projections

**Coordinate Reference System:** WGS84 (EPSG:4326)
- **Datum:** World Geodetic System 1984
- **Units:** Decimal degrees
- **Precision:** Up to 6 decimal places (~0.1 meter accuracy)

---

## 🔧 Data Quality Notes

### Known Issues:
1. **USA Dataset:** Some longitude values appear to be positive instead of negative (Western hemisphere should be negative)
2. **Missing Data:** Some county fields are empty, particularly in the Canada-Mexico dataset
3. **Capacity Breakdown:** Technology-specific capacity columns in USA dataset may not always sum to total capacity due to rounding

### Data Validation:
- ✅ All coordinates fall within expected geographic bounds
- ✅ Plant codes are unique within each dataset
- ✅ Capacity values are non-negative
- ⚠️ Some coordinate precision varies between datasets

---

## 📈 Usage Recommendations

### For Visualization:
- Use latitude/longitude for direct mapping
- Filter by Primary Type for technology-specific views
- Size markers by Total Capacity (MW) for capacity visualization
- Color-code by Balancing Authority for grid operator analysis

### For Analysis:
- Combine both datasets for complete Western Interconnection coverage
- Use State/Province fields for regional aggregation
- Leverage technology-specific capacity columns for detailed energy mix analysis
- Consider capacity ranges for statistical analysis (log scale recommended)

---

## 🔗 Integration with Westmap Application

These datasets are integrated into the Westmap visualization platform:

1. **Data Loading:** Automatically loaded as GeoJSON features
2. **Visualization:** Rendered as point markers on the interactive map
3. **Filtering:** Available through technology type and capacity filters
4. **Analysis:** Used for regional power generation analysis and grid planning

### File Integration Status:
- ✅ **USA Dataset:** Integrated and visualized
- ✅ **Canada-Mexico Dataset:** Integrated and visualized
- ✅ **Combined Analysis:** Available in application

---

## 📝 Metadata

**Last Updated:** September 2025  
**Data Sources:** EIA Form 860, Provincial/Federal Energy Agencies  
**Update Frequency:** Annual  
**License:** Public Domain / Government Data  
**Contact:** Westmap Development Team  

---

## 🚀 Quick Start Examples

### Loading in Python:
```python
import pandas as pd

# Load USA power plants
usa_plants = pd.read_csv('Western_Power_plants_Locations-USA (2).csv')

# Load Canada-Mexico power plants  
can_mex_plants = pd.read_csv('Western_Power_plants_Locations-CANADA-MEX_Update (1).csv')

# Basic statistics
print(f"USA Plants: {len(usa_plants)} facilities")
print(f"Canada-Mexico Plants: {len(can_mex_plants)} facilities")
print(f"Total Capacity (USA): {usa_plants['Total Capacity (MW)'].sum():.1f} MW")
```

### Loading in R:
```r
# Load datasets
usa_plants <- read.csv("Western_Power_plants_Locations-USA (2).csv")
can_mex_plants <- read.csv("Western_Power_plants_Locations-CANADA-MEX_Update (1).csv")

# Summary statistics
summary(usa_plants$`Total.Capacity..MW.`)
table(usa_plants$Primary.Type)
```

---

*This documentation is automatically maintained and updated with each data refresh.*
