import os.path as opath
import csv, pickle
#
from __path_organizer import ef_dpath, pf_dpath

DEPOT, DEPOT_LAT, DEPOT_LNG = 'Keppel Logistics', 1.371509, 103.909272
SEC60 = 60


def get_loctt():
    csv_fpath = opath.join(pf_dpath, 'LocationTravelTime.csv')
    pkl_fpath = opath.join(pf_dpath, 'LocationTravelTime.pkl')
    if opath.exists(pkl_fpath):
        with open(pkl_fpath, 'rb') as fp:
            loctt = pickle.load(fp)
        return loctt
    #
    locations, wholeLocPairs = get_locations_locPairs()
    handled_locPairs = set()
    with open(csv_fpath) as r_csvfile:
        reader = csv.DictReader(r_csvfile)
        for row in reader:
            Location0, Location1 = [row[cn] for cn in ['Location0', 'Location1']]
            handled_locPairs.add((Location0, Location1))
    target_locPairs = wholeLocPairs.difference(handled_locPairs)
    if len(target_locPairs) != 0:
        get_loctt_googleMap()
    loctt = {}
    with open(csv_fpath) as r_csvfile:
        reader = csv.DictReader(r_csvfile)
        for row in reader:
            Location0, Location1 = [row[cn] for cn in ['Location0', 'Location1']]
            Duration = eval(row['Duration'])
            loctt[Location0, Location1] = Duration
    with open(pkl_fpath, 'wb') as fp:
        pickle.dump(loctt, fp)
    #
    return loctt


def get_loctt_googleMap():
    import googlemaps
    googleKeys = [
        'AIzaSyAQYLeLHyJvNVC7uIbHmnvf7x9XC6murmk',
        'AIzaSyDCiqj9QQ-lXWGmzxXM0j-Gbeo_BRlsd0g',
        'AIzaSyCsrxK4ZuxQAsGYt3RNHLeGfEFHwq-GIEU',
        'AIzaSyB2mRWLDgNcAi99A8wGQXqCqecHjihzEa0',
        'AIzaSyBBLZdM06kYG7hu1EmJXGjyL4116Ss1ZBw',
        'AIzaSyCGjnummRXX6WC9CLNWLBJcEhwjSxMvu3U',
        'AIzaSyBDHNYt3PaJ6ctjpSEkcy6zGex2YkVvrpI',
        'AIzaSyDjIRqts0oqVrhPeeO70HaSWx3HhN78s7A',
        'AIzaSyDRwbDVp0pXjTjS838gcM9oy_ssn_QvvdA',
        'AIzaSyDjTCRlNSyWMf2x3INTA-pwADWTwIJFtaE',
        'AIzaSyDdKkkcgB5B2RpqD_GOaNeco4uxUXwodDk'
    ]
    numTrial = 0
    locations, wholeLocPairs = get_locations_locPairs()
    while True:
        print('numTrial', numTrial)
        try:
            googleKey = googleKeys[numTrial % len(googleKeys)]
            gmaps = googlemaps.Client(key=googleKey)
            #
            csv_ofpath = opath.join(pf_dpath, 'LocationTravelTime.csv')
            if not opath.exists(csv_ofpath):
                with open(csv_ofpath, 'w') as w_csvfile:
                    writer = csv.writer(w_csvfile, lineterminator='\n')
                    new_header = ['Lat0', 'Lng0',
                                  'Lat1', 'Lng1',
                                  'Duration',
                                  'Location0', 'District0',
                                  'Location1', 'District1']
                    writer.writerow(new_header)
                target_locPairs = wholeLocPairs
            else:
                handled_locPairs = set()
                with open(csv_ofpath) as r_csvfile:
                    reader = csv.DictReader(r_csvfile)
                    for row in reader:
                        Location0, Location1 = [row[cn] for cn in ['Location0', 'Location1']]
                        handled_locPairs.add((Location0, Location1))
                target_locPairs = wholeLocPairs.difference(handled_locPairs)
            print('num. of remaining pairs', len(target_locPairs))
            for Lat0, Lng0, Location0, District0 in locations:
                for Lat1, Lng1, Location1, District1 in locations:
                    if Location0 == Location1:
                        duration = 0.0
                    else:
                        if (Location0, Location1) not in target_locPairs:
                            continue
                        loc0, loc1 = tuple([Lat0, Lng0]), tuple([Lat1, Lng1])
                        res = gmaps.distance_matrix(loc0, loc1,
                                                    mode="driving")
                        elements = res['rows'][0]['elements']
                        try:
                            duration = elements[0]['duration']['value'] / SEC60
                        except KeyError:
                            print(Location0, Location1)
                            continue
                    new_row = [Lat0, Lng0,
                               Lat1, Lng1,
                               duration,
                               Location0, District0,
                               Location1, District1]
                    with open(csv_ofpath, 'a') as w_csvfile:
                        writer = csv.writer(w_csvfile, lineterminator='\n')
                        writer.writerow(new_row)
        except googlemaps.exceptions.Timeout:
            numTrial += 1
            continue
        break



def get_loctd_googleMap():
    import googlemaps
    googleKeys = [
        'AIzaSyAQYLeLHyJvNVC7uIbHmnvf7x9XC6murmk',
        'AIzaSyDCiqj9QQ-lXWGmzxXM0j-Gbeo_BRlsd0g',
        'AIzaSyCsrxK4ZuxQAsGYt3RNHLeGfEFHwq-GIEU',
        'AIzaSyB2mRWLDgNcAi99A8wGQXqCqecHjihzEa0',
        'AIzaSyBBLZdM06kYG7hu1EmJXGjyL4116Ss1ZBw',
        'AIzaSyCGjnummRXX6WC9CLNWLBJcEhwjSxMvu3U',
        'AIzaSyBDHNYt3PaJ6ctjpSEkcy6zGex2YkVvrpI',
        'AIzaSyDjIRqts0oqVrhPeeO70HaSWx3HhN78s7A',
        'AIzaSyDRwbDVp0pXjTjS838gcM9oy_ssn_QvvdA',
        'AIzaSyDjTCRlNSyWMf2x3INTA-pwADWTwIJFtaE',
        'AIzaSyDdKkkcgB5B2RpqD_GOaNeco4uxUXwodDk'
    ]
    numTrial = 0
    locations, wholeLocPairs = get_locations_locPairs()
    while True:
        print('numTrial', numTrial)
        try:
            googleKey = googleKeys[numTrial % len(googleKeys)]
            gmaps = googlemaps.Client(key=googleKey)
            #
            csv_ofpath = opath.join(pf_dpath, 'LocationTravelDistance.csv')
            if not opath.exists(csv_ofpath):
                with open(csv_ofpath, 'w') as w_csvfile:
                    writer = csv.writer(w_csvfile, lineterminator='\n')
                    new_header = ['Lat0', 'Lng0',
                                  'Lat1', 'Lng1',
                                  'distance',
                                  'Location0', 'District0',
                                  'Location1', 'District1']
                    writer.writerow(new_header)
                target_locPairs = wholeLocPairs
            else:
                handled_locPairs = set()
                with open(csv_ofpath) as r_csvfile:
                    reader = csv.DictReader(r_csvfile)
                    for row in reader:
                        Location0, Location1 = [row[cn] for cn in ['Location0', 'Location1']]
                        handled_locPairs.add((Location0, Location1))
                target_locPairs = wholeLocPairs.difference(handled_locPairs)
            print('num. of remaining pairs', len(target_locPairs))
            for Lat0, Lng0, Location0, District0 in locations:
                for Lat1, Lng1, Location1, District1 in locations:
                    if Location0 == Location1:
                        distance = 0.0
                    else:
                        if (Location0, Location1) not in target_locPairs:
                            continue
                        loc0, loc1 = tuple([Lat0, Lng0]), tuple([Lat1, Lng1])
                        res = gmaps.distance_matrix(loc0, loc1,
                                                    mode="driving")
                        elements = res['rows'][0]['elements']
                        try:
                            distance = elements[0]['distance']['value']
                        except KeyError:
                            print(Location0, Location1)
                            continue
                    new_row = [Lat0, Lng0,
                               Lat1, Lng1,
                               distance,
                               Location0, District0,
                               Location1, District1]
                    with open(csv_ofpath, 'a') as w_csvfile:
                        writer = csv.writer(w_csvfile, lineterminator='\n')
                        writer.writerow(new_row)
        except googlemaps.exceptions.Timeout:
            numTrial += 1
            continue
        break


def get_locations_locPairs():
    csv_fpath = opath.join(ef_dpath, 'LocationPD.csv')
    locations = []
    with open(csv_fpath) as r_csvfile:
        reader = csv.DictReader(r_csvfile)
        for row in reader:
            Lat, Lng = [eval(row[cn]) for cn in ['Lat', 'Lng']]
            Location, District = [row[cn] for cn in ['Location', 'District']]
            locations.append([Lat, Lng, Location, District])
    locations.append([DEPOT_LAT, DEPOT_LNG, DEPOT, None])
    locPairs = set()
    for _, _, loc0, _ in locations:
        for _, _, loc1, _ in locations:
            if loc0 == loc1:
                continue
            locPairs.add((loc0, loc1))
    return locations, locPairs



if __name__ == '__main__':
    print(len(get_loctt()))
