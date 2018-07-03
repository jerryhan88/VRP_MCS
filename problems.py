import os.path as opath
import os
import pickle
#
from sgLocationTravelTime import get_loctt
from sgLocationTravelTime import DEPOT

MAX_NUM_VEH = 50
ALLOWABLE_TT = 5 * 60  # min.
TASK_REWARD = 12  # S$
VEH_COST = 200  # S$


def mrtS1():
    dplym_fpath = opath.join('_temp', 'dplym_%s.pkl' % 'mrtS1-dt80')

    _, problemName = opath.basename(dplym_fpath)[:-len('.pkl')].split('_')
    with open(dplym_fpath, 'rb') as fp:
        _, task_ppdp = pickle.load(fp)
    prmt = convert_prob2prmt(problemName, task_ppdp)
    return prmt


def convert_prob2prmt(problemName, task_ppdp, numVeh=None):
    lid_loc, loc_lid = {}, {}
    lid_loc[0] = DEPOT
    loc_lid[DEPOT] = 0
    numLocs = len(lid_loc)
    for locs in task_ppdp:
        for loc in locs:
            if loc not in loc_lid:
                lid_loc[numLocs] = loc
                loc_lid[loc] = numLocs
                numLocs += 1
    #
    loctt = get_loctt()
    travel_time = [[0.0 for _ in range(numLocs)] for _ in range(numLocs)]
    for pp0, dp0 in task_ppdp:
        travel_time[loc_lid[pp0]][loc_lid[dp0]] = loctt[pp0, dp0]
        travel_time[loc_lid[dp0]][loc_lid[pp0]] = loctt[dp0, pp0]
        #
        travel_time[loc_lid[DEPOT]][loc_lid[pp0]] = loctt[DEPOT, pp0]
        travel_time[loc_lid[pp0]][loc_lid[DEPOT]] = loctt[pp0, DEPOT]
        travel_time[loc_lid[DEPOT]][loc_lid[dp0]] = loctt[DEPOT, dp0]
        travel_time[loc_lid[dp0]][loc_lid[DEPOT]] = loctt[dp0, DEPOT]
        for pp1, dp1 in task_ppdp:
            for loc0, loc1 in [(pp0, pp1), (pp0, dp1),
                               (dp0, pp1), (dp0, dp1)]:
                travel_time[loc_lid[loc0]][loc_lid[loc1]] = loctt[loc0, loc1]
                travel_time[loc_lid[loc1]][loc_lid[loc0]] = loctt[loc1, loc0]
    tasks = []
    for i, (loc0, loc1) in enumerate(task_ppdp):
        tasks.append((i, loc_lid[loc0], loc_lid[loc1]))
    #
    T = list(range(len(tasks)))
    iPs, iMs = list(zip(*[(tasks[i][1], tasks[i][2]) for i in T]))
    P, D = set(), set()
    _N = {}
    for i in T:
        P.add('p%d' % i)
        D.add('d%d' % i)
        #
        _N['p%d' % i] = iPs[i]
        _N['d%d' % i] = iMs[i]
    N = set(_N.keys())
    n0 = 'n0'
    t_ij = {}
    for i in _N:
        t_ij[n0, i] = travel_time[0][_N[i]]
        t_ij[i, n0] = travel_time[_N[i]][0]
        for j in _N:
            t_ij[i, j] = travel_time[_N[i]][_N[j]]
    V = list(range(numVeh if numVeh else MAX_NUM_VEH))
    _delta = ALLOWABLE_TT
    cP, cC = TASK_REWARD, VEH_COST
    #
    return {'problemName': problemName,
            'V': V,
            'T': T, 'P': P, 'D': D, 'N': N,
            'n0': n0,
            't_ij': t_ij, '_delta': _delta,
            'cP': cP, 'cC': cC
    }


if __name__ == '__main__':
    mrtS1()
