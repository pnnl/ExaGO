# Powerlines WUS CAN SGCA Dataset Documentation

## Overview
This dataset contains powerline corridors in the Western United States and Canada, compiled from 22 different source data layers. The data was originally created for the Conservation Assessment of Greater Sage-grouse and Sagebrush Habitats but can be used for any purpose requiring general powerline location information.

## Dataset Information
- **Title**: Powerline Corridors in the Western United States and Canada
- **Author**: Sean P. Finn, USGS Forest and Rangeland Ecosystem Science Center
- **Publication Date**: March 1, 2023 (originally compiled 2004)
- **DOI**: https://doi.org/10.5066/P915FITQ
- **Total Features**: 32,117 powerline segments

## Geographic Coverage
- **Region**: Western United States and Canada
- **Bounding Box**: 
  - West: -136.8556°
  - East: -88.2494°
  - North: 61.7607°
  - South: 30.0784°
- **States/Provinces Covered**: Washington, Oregon, California, Idaho, Nevada, Utah, Arizona, Montana, Wyoming, Colorado, New Mexico, North Dakota, South Dakota, British Columbia, Alberta, Saskatchewan

## File Structure

### Shapefile Components
The dataset consists of several files that together form an ESRI Shapefile:

#### Core Shapefile Files:
1. **`powerlines_WUS_CAN_sgca.shp`** - Main geometry file containing the actual powerline vector data
2. **`powerlines_WUS_CAN_sgca.shx`** - Shape index file for spatial indexing
3. **`powerlines_WUS_CAN_sgca.dbf`** - Attribute database file containing feature attributes
4. **`powerlines_WUS_CAN_sgca.prj`** - Projection file defining coordinate reference system
5. **`powerlines_WUS_CAN_sgca.sbn`** - Spatial index file for ArcGIS (binary)
6. **`powerlines_WUS_CAN_sgca.sbx`** - Spatial index file for ArcGIS (index)

#### Metadata and Documentation:
7. **`powerlines_WUS_CAN_sgca.shp.xml`** - Comprehensive XML metadata file with detailed dataset information

### Converted Files
8. **`powerlines_WUS_CAN_sgca.csv`** - Basic CSV export of attribute data
9. **`powerlines_WUS_CAN_sgca_wgs84.csv`** - Enhanced CSV with WGS84 coordinates and additional fields for visualization

## Coordinate Reference System
- **Original Projection**: NAD 1927 Albers Conical Equal Area
- **Parameters**:
  - Standard Parallels: 38.0°, 41.0°
  - Central Meridian: -114.0°
  - Latitude of Origin: 34.0°
  - Units: Meters
- **Converted Version**: WGS84 (EPSG:4326) for web mapping compatibility

## Data Attributes

### Original Shapefile Attributes:
| Field | Type | Description |
|-------|------|-------------|
| `FID` | Integer | Unique feature identifier |
| `Shape` | Geometry | Polyline geometry of powerlines |
| `LENGTH` | Numeric | Length of powerline segment (units unknown) |
| `CA_POWER` | Text(50) | Type of linear feature (POWER, PIPE, OTHER, or blank) |
| `SOURCEFILE` | Text(50) | Original source dataset identifier |

### Enhanced CSV Attributes (WGS84 version):
| Field | Type | Description |
|-------|------|-------------|
| `wkt` | Text | Well-Known Text geometry representation |
| `flow capacity` | Numeric | Power flow capacity (default: 100) |
| `pf`, `qf`, `pt`, `qt` | Numeric | Power flow parameters (from/to active/reactive power) |
| `kilovolt` | Numeric | Voltage level (default: 230 kV) |
| `line_name` | Text | Generated line identifier |
| `source` | Text | Source type (typically "POWER") |
| `target` | Text | Target/destination identifier |
| `actual flow` | Numeric | Current power flow (default: 0) |
| `data_source` | Text | Original dataset identifier |

## Data Quality and Limitations

### Accuracy
- **Positional Accuracy**: Varies by source data; recommended for use at scales of 1:100,000 and smaller
- **Temporal Currency**: Data compiled as of May 15, 2004
- **Completeness**: NOT intended to be a complete representation of all powerlines on the ground
- **Focus**: Primarily depicts higher voltage, long-distance transmission lines; some segments may represent lower voltage distribution lines

### Source Data Reliability
The dataset was compiled from 22 different sources with varying quality and accuracy levels. Users should consider this variability when using the data for analysis.

## Data Sources
The dataset combines powerline data from multiple sources including:
- National Park Service facilities
- Bureau of Land Management field offices  
- State utility agencies
- US Census Bureau TIGER data
- Digital Chart of the World
- Utility companies (PacifiCorp, etc.)
- University research projects

For a complete list of 22 source datasets, refer to the XML metadata file.

## Usage in Visualization
This dataset is used in the ExaGO visualization platform to display transmission line infrastructure on the map. The WGS84 CSV version provides compatibility with web mapping frameworks like Deck.gl and React-map-gl.

### Visualization Features:
- Power flow visualization along transmission corridors
- Interactive selection and filtering by voltage levels
- Integration with generation and load data
- Geographic context for power system analysis

## Data Processing
The original shapefile data has been processed to create web-compatible formats:
1. Reprojection from NAD 1927 Albers to WGS84
2. Conversion to CSV format with WKT geometry
3. Addition of power system attributes for visualization
4. Integration with ExaGO power system analysis results

## Disclaimer
This data is provided "as is" without warranty. Users should read the complete metadata and understand data limitations before use. The dataset is public domain and freely redistributable with proper attribution to USGS.

## Citation
Finn, Sean P., 2023, Powerline Corridors in the Western United States and Canada: U.S. Geological Survey data release, https://doi.org/10.5066/P915FITQ

## Contact
For questions about this dataset:
- **Organization**: USGS Forest and Rangeland Ecosystem Science Center
- **Email**: fresc_outreach@usgs.gov
- **Phone**: 541-750-1030
- **Address**: 777 NW 9th St, Suite 400, Corvallis, OR 97330
