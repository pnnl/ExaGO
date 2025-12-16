import pandas as pd
import json
from shapely.geometry import shape
import csv
import sys

import geopandas as gpd
from shapely.geometry import Point
import requests, zipfile, io, os, tempfile

URL = "https://www2.census.gov/geo/tiger/TIGER2024/COUNTY/tl_2024_us_county.zip"


def load_us_counties(force_download=False):
    cache_dir = os.path.join(tempfile.gettempdir(), "us_counties")
    shp_path = os.path.join(cache_dir, "tl_2024_us_county.shp")
    if force_download or not os.path.exists(shp_path):
        os.makedirs(cache_dir, exist_ok=True)
        r = requests.get(URL, timeout=60)
        r.raise_for_status()
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall(cache_dir)
    gdf = gpd.read_file(shp_path)
    gdf = gdf.to_crs(epsg=4326)
    print(gdf.columns)
    return gdf[["NAME", "STATEFP", "GEOID", "geometry"]]


counties = load_us_counties()
STATE_FIPS = {
    "01": ("AL", "Alabama"),
    "02": ("AK", "Alaska"),
    "04": ("AZ", "Arizona"),
    "05": ("AR", "Arkansas"),
    "06": ("CA", "California"),
    "08": ("CO", "Colorado"),
    "09": ("CT", "Connecticut"),
    "10": ("DE", "Delaware"),
    "11": ("DC", "District of Columbia"),
    "12": ("FL", "Florida"),
    "13": ("GA", "Georgia"),
    "15": ("HI", "Hawaii"),
    "16": ("ID", "Idaho"),
    "17": ("IL", "Illinois"),
    "18": ("IN", "Indiana"),
    "19": ("IA", "Iowa"),
    "20": ("KS", "Kansas"),
    "21": ("KY", "Kentucky"),
    "22": ("LA", "Louisiana"),
    "23": ("ME", "Maine"),
    "24": ("MD", "Maryland"),
    "25": ("MA", "Massachusetts"),
    "26": ("MI", "Michigan"),
    "27": ("MN", "Minnesota"),
    "28": ("MS", "Mississippi"),
    "29": ("MO", "Missouri"),
    "30": ("MT", "Montana"),
    "31": ("NE", "Nebraska"),
    "32": ("NV", "Nevada"),
    "33": ("NH", "New Hampshire"),
    "34": ("NJ", "New Jersey"),
    "35": ("NM", "New Mexico"),
    "36": ("NY", "New York"),
    "37": ("NC", "North Carolina"),
    "38": ("ND", "North Dakota"),
    "39": ("OH", "Ohio"),
    "40": ("OK", "Oklahoma"),
    "41": ("OR", "Oregon"),
    "42": ("PA", "Pennsylvania"),
    "44": ("RI", "Rhode Island"),
    "45": ("SC", "South Carolina"),
    "46": ("SD", "South Dakota"),
    "47": ("TN", "Tennessee"),
    "48": ("TX", "Texas"),
    "49": ("UT", "Utah"),
    "50": ("VT", "Vermont"),
    "51": ("VA", "Virginia"),
    "53": ("WA", "Washington"),
    "54": ("WV", "West Virginia"),
    "55": ("WI", "Wisconsin"),
    "56": ("WY", "Wyoming"),
    "72": ("PR", "Puerto Rico"),
}


def county_state_from_latlon(lat, lon, counties_gdf=counties):
    pt = Point(lon, lat)
    # quick bounding-box filter
    candidates = list(counties_gdf.sindex.intersection(pt.bounds))
    if not candidates:
        print("latlon not found:", lat, lon)
        return None, None
    subset = counties_gdf.iloc[candidates]
    match = subset[subset.contains(pt)]

    if match.empty:
        print("latlon not found:", lat, lon)
        return None, None
    row = match.iloc[0]
    state = STATE_FIPS.get(row["STATEFP"], ("", ""))[1]
    county = row["NAME"]
    return county, state


def getGeneration(data):
    Geni = None
    Gens = []
    minPg = 1000.0
    maxPg = 0.0
    Pgcoal = 0.0
    Pghydro = 0.0
    Pgnuclear = 0.0
    Pgng = 0.0
    Pgsolar = 0.0
    Pgwind = 0.0
    Pgother = 0.0
    Pgcoalcap = 0.0
    Pghydrocap = 0.0
    Pgnuclearcap = 0.0
    Pgngcap = 0.0
    Pgsolarcap = 0.0
    Pgwindcap = 0.0
    Pgothercap = 0.0

    for feature in data["geojsondata"]["features"]:
        if feature["geometry"]["type"] == "Point":
            subst = feature["properties"]
            name = subst["NAME"]
            nbus = subst["nbus"]
            Pg = 0.0
            Pcap = 0.0
            gen_fuel = ""
            ngen = 0
            KV = []
            for bus in subst["bus"]:
                KV.append(bus["BASE_KV"])

                for gen in bus["gen"]:
                    Pg += gen["GEN_STATUS"] * gen["PG"]
                    Pcap += gen["GEN_STATUS"] * gen["PMAX"]
                    gen_fuel = gen["GEN_FUEL"].lower()

                    if gen_fuel == "wind":
                        Pgwind += gen["GEN_STATUS"] * gen["PG"]
                        Pgwindcap += gen["GEN_STATUS"] * gen["PMAX"]
                    elif gen_fuel == "solar":
                        Pgsolar += gen["GEN_STATUS"] * gen["PG"]
                        Pgsolarcap += gen["GEN_STATUS"] * gen["PMAX"]
                    elif gen_fuel == "coal":
                        Pgcoal += gen["GEN_STATUS"] * gen["PG"]
                        Pgcoalcap += gen["GEN_STATUS"] * gen["PMAX"]
                    elif gen_fuel == "nuclear":
                        Pgnuclear += gen["GEN_STATUS"] * gen["PG"]
                        Pgnuclearcap += gen["GEN_STATUS"] * gen["PMAX"]
                    elif gen_fuel == "hydro":
                        Pghydro += gen["GEN_STATUS"] * gen["PG"]
                        Pghydrocap += gen["GEN_STATUS"] * gen["PMAX"]
                    elif gen_fuel == "ng":
                        Pgng += gen["GEN_STATUS"] * gen["PG"]
                        Pgngcap += gen["GEN_STATUS"] * gen["PMAX"]
                    else:
                        Pgother += gen["GEN_STATUS"] * gen["PG"]
                        Pgothercap += gen["GEN_STATUS"] * gen["PMAX"]
                    ngen += 1

            if ngen:
                color = ""
                if gen_fuel == "wind":
                    color = "green"
                elif gen_fuel == "solar":
                    color = "yellow"
                elif gen_fuel == "coal":
                    color = "gray"
                elif gen_fuel == "nuclear":
                    color = "red"
                elif gen_fuel == "hydro":
                    color = "blue"
                elif gen_fuel == "ng":
                    color = "orange"
                else:
                    color = "black"

                if Pg <= minPg:
                    minPg = Pg
                if Pg >= maxPg:
                    maxPg = Pg
                geo = shape(feature["geometry"])

                coords = geo.coords[0]  # Get first coordinate pair
                lat, lon = coords[1], coords[0]  # GeoJSON uses (lon, lat) order
                county_info = county_state_from_latlon(lat, lon)

                Geni = {
                    "coordinates": geo.wkt,
                    "Power generated": Pg,
                    "Power capacity": Pcap,
                    "Remaining Capacity": abs(Pcap - Pg),
                    "KVlevels": set(KV),
                    "color": color,
                    "generation name": name,
                    "number of buses": nbus,
                    "generation type": gen_fuel,
                    "county": county_info[0],
                    "state": county_info[1],
                }
                Gens.append(Geni)

    keys = Gens[0].keys()
    with open("generation.csv", "w", newline="") as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(Gens)


def getBus(data):
    # f = open('case_ACTIVSg10k point.json')
    # data = json.load(f)
    points = []
    for feature in data["geojsondata"]["features"]:
        if feature["geometry"]["type"] == "Point":
            geo = shape(feature["geometry"])
            p = feature["properties"]
            # kvlevels = p["KVlevels"].replace("[", "{")
            # kvlevels = kvlevels.replace("]", "}")

            coords = geo.coords[0]  # Get first coordinate pair
            lat, lon = coords[1], coords[0]  # GeoJSON uses (lon, lat) order
            county_info = county_state_from_latlon(lat, lon)

            points.append(
                {
                    "coordinates": geo.wkt,
                    "bus_name": p["NAME"],
                    "kilovolt levels": set(p["KVlevels"]),
                    "number of buses": p["nbus"],
                    "vm": p["Vm"],
                    "county": county_info[0],
                    "state": county_info[1],
                }
            )

    keys = points[0].keys()
    with open("bus.csv", "w", newline="") as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(points)


def getLine(data):
    # f = open('case_ACTIVSg10k line.json')
    # data = json.load(f)
    lines = []
    for feature in data["geojsondata"]["features"]:
        if feature["geometry"]["type"] == "LineString":
            geo = shape(feature["geometry"])
            p = feature["properties"]
            # kvlevels = p["KVlevels"].replace("[", "{")
            # kvlevels = kvlevels.replace("]", "}")
            x = p["NAME"].split(" -- ")

            if p["PF"] > 0:
                lines.append(
                    {
                        "wkt": geo.wkt,
                        "flow capacity": p["RATE_A"],
                        "pf": p["PF"],
                        "qf": p["QF"],
                        "pt": p["PT"],
                        "qt": p["QT"],
                        "kilovolt": p["KV"],
                        "line_name": p["NAME"],
                        "srouce": x[0],
                        "target": x[1],
                        "actual flow": abs(p["PF"]),
                        "loading percent": abs(p["PF"]) / p["RATE_A"] * 100,
                    }
                )
            else:
                lines.append(
                    {
                        "wkt": geo.wkt,
                        "flow capacity": p["RATE_A"],
                        "pf": p["PF"],
                        "qf": p["QF"],
                        "pt": p["PT"],
                        "qt": p["QT"],
                        "kilovolt": p["KV"],
                        "line_name": p["NAME"],
                        "srouce": x[1],
                        "target": x[0],
                        "actual flow": abs(p["PF"]),
                        "loading percent": abs(p["PF"]) / p["RATE_A"] * 100,
                    }
                )

    keys = lines[0].keys()

    with open("transmission_line.csv", "w", newline="") as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(lines)


def main():
    filename = sys.argv[1]
    f = open(filename)
    data = json.load(f)
    getGeneration(data)
    getBus(data)
    getLine(data)


if __name__ == "__main__":
    main()
