import pandas as pd
import json
from shapely.geometry import shape
import csv
import sys


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

    for feature in data['geojsondata']['features']:
        if feature['geometry']['type'] == 'Point':
            subst = feature['properties']
            name = subst['NAME']
            nbus = subst['nbus']
            Pg = 0.0
            Pcap = 0.0
            gen_fuel = ''
            ngen = 0
            KV = []
            for bus in subst['bus']:
                KV.append(bus['BASE_KV'])

                for gen in bus['gen']:
                    Pg += gen['GEN_STATUS'] * gen['PG']
                    Pcap += gen['GEN_STATUS'] * gen['PMAX']
                    gen_fuel = gen['GEN_FUEL'].lower()

                    if gen_fuel == 'wind':
                        Pgwind += gen['GEN_STATUS'] * gen['PG']
                        Pgwindcap += gen['GEN_STATUS'] * gen['PMAX']
                    elif gen_fuel == 'solar':
                        Pgsolar += gen['GEN_STATUS'] * gen['PG']
                        Pgsolarcap += gen['GEN_STATUS'] * gen['PMAX']
                    elif gen_fuel == 'coal':
                        Pgcoal += gen['GEN_STATUS'] * gen['PG']
                        Pgcoalcap += gen['GEN_STATUS'] * gen['PMAX']
                    elif gen_fuel == 'nuclear':
                        Pgnuclear += gen['GEN_STATUS'] * gen['PG']
                        Pgnuclearcap += gen['GEN_STATUS'] * gen['PMAX']
                    elif gen_fuel == 'hydro':
                        Pghydro += gen['GEN_STATUS'] * gen['PG']
                        Pghydrocap += gen['GEN_STATUS'] * gen['PMAX']
                    elif gen_fuel == 'ng':
                        Pgng += gen['GEN_STATUS'] * gen['PG']
                        Pgngcap += gen['GEN_STATUS'] * gen['PMAX']
                    else:
                        Pgother += gen['GEN_STATUS'] * gen['PG']
                        Pgothercap += gen['GEN_STATUS'] * gen['PMAX']
                    ngen += 1

            if ngen:
                color = ''
                if gen_fuel == 'wind':
                    color = 'green'
                elif gen_fuel == 'solar':
                    color = 'yellow'
                elif gen_fuel == 'coal':
                    color = 'gray'
                elif gen_fuel == 'nuclear':
                    color = 'red'
                elif gen_fuel == 'hydro':
                    color = 'blue'
                elif gen_fuel == 'ng':
                    color = 'orange'
                else:
                    color = 'black'

                if Pg <= minPg:
                    minPg = Pg
                if Pg >= maxPg:
                    maxPg = Pg
                geo = shape(feature["geometry"])
                Geni = {"coordinates": geo.wkt, "Power generated": Pg, "Power capacity": Pcap, "KVlevels": set(
                    KV), "color": color, "generation name": name, "number of buses": nbus, "generation type": gen_fuel}
                Gens.append(Geni)

    keys = Gens[0].keys()
    with open('generation.csv', 'w', newline='') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(Gens)


def getBus(data):
    # f = open('case_ACTIVSg10k point.json')
    # data = json.load(f)
    points = []
    for feature in data['geojsondata']['features']:
        if feature['geometry']['type'] == 'Point':
            geo = shape(feature["geometry"])
            p = feature["properties"]
            # kvlevels = p["KVlevels"].replace("[", "{")
            # kvlevels = kvlevels.replace("]", "}")
            points.append({
                "wkt": geo.wkt,
                "bus_name": p["NAME"],
                "kilovolt levels": set(p["KVlevels"]),
                "number of buses": p["nbus"],
                "vm": p["Vm"],
            })

    keys = points[0].keys()
    with open('bus.csv', 'w', newline='') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(points)


def getLine(data):
    # f = open('case_ACTIVSg10k line.json')
    # data = json.load(f)
    lines = []
    for feature in data['geojsondata']['features']:
        if feature['geometry']['type'] == "LineString":
            geo = shape(feature["geometry"])
            p = feature["properties"]
            # kvlevels = p["KVlevels"].replace("[", "{")
            # kvlevels = kvlevels.replace("]", "}")
            x = p["NAME"].split(' -- ')
            if (p['PF'] > 0):
                lines.append({
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


                })
            else:
                lines.append({
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
                })

    keys = lines[0].keys()

    with open('transmission_line.csv', 'w', newline='') as output_file:
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
