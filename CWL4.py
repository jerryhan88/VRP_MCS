import os.path as opath
import multiprocessing
import time
import pickle, csv
import datetime
import numpy as np
from gurobipy import *
#
from _util_cython import gen_cFile

prefix = 'PD_IH'
gen_cFile(prefix)
from PD_IH import run as PD_IH_run
from PD_IH import calc_travelTime


NUM_CORES = multiprocessing.cpu_count()
EPSILON = 0.00001


def write_log(fpath, contents):
    with open(fpath, 'a') as f:
        logContents = '\n\n'
        logContents += '======================================================================================\n'
        logContents += '%s\n' % str(datetime.datetime.now())
        logContents += '%s\n' % contents
        logContents += '======================================================================================\n'
        f.write(logContents)


def itr2file(fpath, contents=[]):
    if not contents:
        if opath.exists(fpath):
            os.remove(fpath)
        with open(fpath, 'wt') as w_csvfile:
            writer = csv.writer(w_csvfile, lineterminator='\n')
            header = ['itrNum',
                      'eliCpuTime', 'eliWallTime',
                      'numCols', 'numTB',
                      'relObjV', 'selBC', 'selBC_RC', 'new_RC_BC']
            writer.writerow(header)
    else:
        with open(fpath, 'a') as w_csvfile:
            writer = csv.writer(w_csvfile, lineterminator='\n')
            writer.writerow(contents)


def res2file(fpath, objV, gap, eliCpuTime, eliWallTime):
    with open(fpath, 'wt') as w_csvfile:
        writer = csv.writer(w_csvfile, lineterminator='\n')
        header = ['objV', 'Gap', 'eliCpuTime', 'eliWallTime']
        writer.writerow(header)
        writer.writerow([objV, gap, eliCpuTime, eliWallTime])


def generate_RMP(prmt, add_inputs):
    C, p_c, e_ci = list(map(add_inputs.get, ['C', 'p_c', 'e_ci']))
    #
    V, T, cC = list(map(prmt.get, ['V', 'T', 'cC']))
    nV = len(V)
    #
    # Define decision variables
    #
    RMP = Model('RMP')
    q_c = {}
    for c in range(len(C)):
        q_c[c] = RMP.addVar(vtype=GRB.BINARY, name="q[%d]" % c)
    RMP.update()
    #
    # Define objective
    #
    obj = LinExpr()
    for c in range(len(C)):
        obj += p_c[c] * q_c[c]
    # for c in range(len(C)):
    #     obj -= cC * q_c[c]
    RMP.setObjective(obj, GRB.MAXIMIZE)
    #
    # Define constrains
    #
    taskAC = {}
    for i in T:  # eq:taskA
        taskAC[i] = RMP.addConstr(quicksum(e_ci[c][i] * q_c[c] for c in range(len(C))) <= 1,
                                  name="taskAC[%d]" % i)
    numVC = RMP.addConstr(quicksum(q_c[c] for c in range(len(C))) <= nV,
                              name="numVC")
    RMP.update()
    #
    return RMP, q_c, taskAC, numVC


def LS_run(prmt, cwl_inputs):
    T, cP, cC, _delta, t_ij = [prmt.get(k) for k in ['T', 'cP', 'cC', '_delta', 't_ij']]
    #
    pi_i, mu = [cwl_inputs.get(k) for k in ['pi_i', 'mu']]
    C, sC, c0 = [cwl_inputs.get(k) for k in ['C', 'sC', 'c0']]
    s_c = cwl_inputs['s_c']
    mT = cwl_inputs['mT']
    #
    Ts0 = C[c0]
    tt_rc_Ts1_seq = []
    for i0 in T:
        if i0 in Ts0:
            continue
        Ts1 = Ts0[:] + [i0]
        if frozenset(tuple(Ts1)) in sC:
            continue
        seq0 = s_c[c0]
        travelTime, seq1 = PD_IH_run(prmt, {'seq0': seq0, 'i0': i0})
        woDetSeq = seq1[1:-1]
        woDetTT = calc_travelTime(woDetSeq, t_ij)
        if woDetTT <= len(woDetSeq) * mT * 0.25 and travelTime <= _delta:
        # if travelTime <= _delta:
            vec = [0 for _ in range(len(T))]
            for i in Ts1:
                vec[i] = 1
            rc = cP * len(Ts1) -cC - (np.array(vec ) * np.array(pi_i)).sum() - mu
            tt_rc_Ts1_seq.append([travelTime, rc, Ts1, seq1])
    tt_rc_Ts1_seq.sort()
    #
    if 2 <= len(tt_rc_Ts1_seq):
        rc_Ts1_seq = [[rc, Ts1, seq1] for _, rc, Ts1, seq1 in tt_rc_Ts1_seq[:2]]
    else:
        rc_Ts1_seq = [[rc, Ts1, seq1] for _, rc, Ts1, seq1 in tt_rc_Ts1_seq]
    return rc_Ts1_seq


def run(prmt, etc=None):
    startCpuTime, startWallTime = time.clock(), time.time()
    if 'TimeLimit' not in etc:
        etc['TimeLimit'] = 1e400
    etc['startTS'] = startCpuTime
    etc['startCpuTime'] = startCpuTime
    etc['startWallTime'] = startWallTime
    itr2file(etc['itrFileCSV'])
    #
    cwl_inputs = {}
    T, cP, cC, n0, t_ij = [prmt.get(k) for k in ['T', 'cP', 'cC', 'n0', 't_ij']]
    #
    C, sC, p_c, e_ci = [], set(), [], []
    TB = set()
    s_c = {}
    for i in T:
        c = len(C)
        iP, iM = 'p%d' % i, 'd%d' % i
        s_c[c] = [n0, iP, iM, n0]
        Ts = [i]
        C.append(Ts)
        sC.add(frozenset(tuple(Ts)))
        #
        p_c.append(cP - cC)
        # p_c.append(cP)
        #
        vec = [0 for _ in range(len(T))]
        vec[i] = 1
        e_ci.append(vec)
    #
    cwl_inputs['C'] = C
    cwl_inputs['sC'] = sC
    cwl_inputs['p_c'] = p_c
    cwl_inputs['e_ci'] = e_ci
    cwl_inputs['TB'] = TB
    cwl_inputs['s_c'] = s_c
    cwl_inputs['mT'] = max(list(t_ij.values()))
    #
    RMP, q_c, taskAC, numVC = generate_RMP(prmt, cwl_inputs)
    #
    counter, is_terminated = 0, False
    while True:
        if len(C) == len(T) ** 2 - 1:
            break
        LRMP = RMP.relax()
        LRMP.setParam('Threads', NUM_CORES)
        LRMP.setParam('OutputFlag', False)
        LRMP.optimize()
        if LRMP.status == GRB.Status.INFEASIBLE:
            logContents = 'Relaxed model is infeasible!!\n'
            logContents += 'No solution!\n'
            write_log(etc['logFile'], logContents)
            #
            LRMP.write('%s.lp' % prmt['problemName'])
            LRMP.computeIIS()
            LRMP.write('%s.ilp' % prmt['problemName'])
            assert False
        #
        pi_i = [LRMP.getConstrByName("taskAC[%d]" % i).Pi for i in T]
        mu = LRMP.getConstrByName("numVC").Pi
        cwl_inputs['pi_i'] = pi_i
        cwl_inputs['mu'] = mu
        #
        c0, minRC = -1, 1e400
        for rc, c in [(LRMP.getVarByName("q[%d]" % c).RC, c) for c in range(len(C))]:
            if c in TB:
                continue
            if rc < -cC:
                TB.add(c)
                continue
            if rc < minRC:
                minRC = rc
                c0 = c
        if c0 == -1:
            break
        cwl_inputs['c0'] = c0
        #
        rc_Ts1_seq = LS_run(prmt, cwl_inputs)
        rc_Ts1 = [[o[0], o[1]] for o in rc_Ts1_seq]
        if time.clock() - etc['startTS'] > etc['TimeLimit']:
            break
        #
        eliCpuTimeP, eliWallTimeP = time.clock() - etc['startCpuTime'], time.time() - etc['startWallTime']
        itr2file(etc['itrFileCSV'], [counter, '%.2f' % eliCpuTimeP, '%.2f' % eliWallTimeP,
                                     len(cwl_inputs['C']), len(cwl_inputs['TB']),
                                     '%.2f' % LRMP.objVal, C[c0], '%.2f' % minRC, str(rc_Ts1)])
        # if len(rc_Ts1_seq) == 0:
        #     TB.add(c0)
        # else:
        #     # is_updated = False
        for rc, Ts1, seq in rc_Ts1_seq:
            if rc < -cC:
                continue
            is_updated = True
            vec = [0 for _ in range(len(T))]
            for i in Ts1:
                vec[i] = 1
            p = rc + (np.array(vec) * np.array(pi_i)).sum() + mu
            C, p_c, e_ci, sC = list(map(cwl_inputs.get, ['C', 'p_c', 'e_ci', 'sC']))
            e_ci.append(vec)
            p_c.append(p)
            #
            col = Column()
            for i in range(len(T)):
                if e_ci[len(C)][i] > 0:
                    col.addTerms(e_ci[len(C)][i], taskAC[i])
            col.addTerms(1, numVC)
            #
            q_c[len(C)] = RMP.addVar(obj=p_c[len(C)], vtype=GRB.BINARY, name="q[%d]" % len(C), column=col)
            s_c[len(C)] = seq
            C.append(Ts1)
            sC.add(frozenset(tuple(Ts1)))
            RMP.update()
            #
        # if not is_updated:
        TB.add(c0)
        if len(C) == len(TB):
            break
        counter += 1
    #
    # Handle termination
    #
    RMP.setParam('Threads', NUM_CORES)
    RMP.optimize()
    RMP.write('temp.lp')
    #
    if etc and RMP.status != GRB.Status.INFEASIBLE:
        assert 'solFilePKL' in etc
        assert 'solFileCSV' in etc
        assert 'solFileTXT' in etc
        #
        q_c = [RMP.getVarByName("q[%d]" % c).x for c in range(len(C))]
        chosenC = [(C[c], '%.2f' % q_c[c]) for c in range(len(C)) if q_c[c] > 0.5]
        with open(etc['solFileTXT'], 'w') as f:
            endCpuTime, endWallTime = time.clock(), time.time()
            eliCpuTime = endCpuTime - etc['startCpuTime']
            eliWallTime = endWallTime - etc['startWallTime']
            logContents = 'Summary\n'
            logContents += '\t Cpu Time: %f\n' % eliCpuTime
            logContents += '\t Wall Time: %f\n' % eliWallTime
            logContents += '\t ObjV: %.3f\n' % RMP.objVal
            logContents += '\t Gap: %.3f\n' % RMP.MIPGap
            logContents += 'Chosen bundles\n'
            logContents += '%s\n' % str(chosenC)
            f.write(logContents)
            f.write('\n')
        #
        res2file(etc['solFileCSV'], RMP.objVal, RMP.MIPGap, eliCpuTime, eliWallTime)
        #
        _q_c = {c: RMP.getVarByName("q[%d]" % c).x for c in range(len(C))}
        sol = {
            'C': C, 'p_c': p_c, 'e_ci': e_ci,
            #
            'q_c': _q_c}
        with open(etc['solFilePKL'], 'wb') as fp:
            pickle.dump(sol, fp)


if __name__ == '__main__':
    from problems import mrtS1
    #
    prmt = mrtS1()
    problemName = prmt['problemName']
    #
    etc = {'solFilePKL': opath.join('_temp', 'sol_%s_CWL4.pkl' % problemName),
           'solFileCSV': opath.join('_temp', 'sol_%s_CWL4.csv' % problemName),
           'solFileTXT': opath.join('_temp', 'sol_%s_CWL4.txt' % problemName),
           'logFile': opath.join('_temp', '%s_CWL4.log' % problemName),
           'itrFileCSV': opath.join('_temp', '%s_itrCWL4.csv' % problemName),
           }
    #
    run(prmt, etc)
