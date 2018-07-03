import os.path as opath
import os
import csv, pickle
import shutil
import pandas as pd
import numpy as np
from functools import reduce
from itertools import chain
from random import normalvariate, seed, sample
#
from __path_organizer import exp_dpath
from problems import convert_prob2prmt
#
from _util_cython import gen_cFile


def gen_prmts(machine_dpath, numVehs=None):
    assert opath.exists(machine_dpath)
    dplym_dpath = opath.join(machine_dpath, 'dplym')
    assert opath.exists(dplym_dpath)
    prmt_dpath = opath.join(machine_dpath, 'prmt')
    if opath.exists(prmt_dpath):
        shutil.rmtree(prmt_dpath)
    os.mkdir(prmt_dpath)
    for fn in os.listdir(dplym_dpath):
        if not fn.endswith('.pkl'):
            continue
        _, problemName = fn[:-len('.pkl')].split('_')
        with open(opath.join(dplym_dpath, fn), 'rb') as fp:
            _, task_ppdp = pickle.load(fp)
        if not numVehs:
            prmt = convert_prob2prmt(problemName, task_ppdp)
            with open(opath.join(prmt_dpath, 'prmt_%s.pkl' % problemName), 'wb') as fp:
                pickle.dump(prmt, fp)
        else:
            assert type(numVehs) == list
            for numVeh in numVehs:
                prmt = convert_prob2prmt(problemName, task_ppdp, numVeh)
                with open(opath.join(prmt_dpath, 'prmt_nv%d_%s.pkl' % (numVeh, problemName)), 'wb') as fp:
                    pickle.dump(prmt, fp)


def run_experiments(machine_num, fixVN=0):
    machine_dpath = opath.join(exp_dpath, 'm%d' % machine_num)
    if not fixVN:
        gen_cFile('CWL4')
        from CWL4 import run as CWL4_run
        assert opath.exists(machine_dpath)
        gen_prmts(machine_dpath)
    else:
        gen_cFile('CWL4_fixVN')
        from CWL4_fixVN import run as CWL4_run
    prmt_dpath = opath.join(machine_dpath, 'prmt')
    log_dpath = opath.join(machine_dpath, 'log')
    sol_dpath = opath.join(machine_dpath, 'sol')
    for path in [log_dpath, sol_dpath]:
        if opath.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
    problems_fpaths = [opath.join(prmt_dpath, fn) for fn in os.listdir(prmt_dpath)
                       if fn.endswith('.pkl')]
    problems_fpaths.sort()
    for fpath in problems_fpaths:
        with open(fpath, 'rb') as fp:
            prmt = pickle.load(fp)
        problemName = prmt['problemName']
        cwl_no = 4
        etc = {'solFilePKL': opath.join(sol_dpath, 'sol_%s_CWL%d.pkl' % (problemName, cwl_no)),
               'solFileCSV': opath.join(sol_dpath, 'sol_%s_CWL%d.csv' % (problemName, cwl_no)),
               'solFileTXT': opath.join(sol_dpath, 'sol_%s_CWL%d.txt' % (problemName, cwl_no)),
               'logFile': opath.join(log_dpath, '%s_CWL%d.log' % (problemName, cwl_no)),
               'itrFileCSV': opath.join(log_dpath, '%s_itrCWL%d.csv' % (problemName, cwl_no)),
               }
        CWL4_run(prmt, etc)


def summaryNV():
    summaryNV_dpath = opath.join(exp_dpath, '_summaryNV_temp')
    sum_fpath = reduce(opath.join, [summaryNV_dpath, 'expSumNV.csv'])
    with open(sum_fpath, 'w') as w_csvfile:
        writer = csv.writer(w_csvfile, lineterminator='\n')
        header = ['pn', 'numTasks', 'allowTT', 'cpuT', 'objV', 'NV', 'HNT', 'avgTPV']
        writer.writerow(header)
    prmt_dpath = reduce(opath.join, [summaryNV_dpath, 'prmt'])
    sol_dpath = reduce(opath.join, [summaryNV_dpath, 'sol'])
    log_dpath = reduce(opath.join, [summaryNV_dpath, 'log'])
    aprc = 'CWL4'
    for fn in os.listdir(prmt_dpath):
        if not fn.endswith('.pkl'): continue
        _, prefix = fn[:-len('.pkl')].split('_')
        prmt_fpath = opath.join(prmt_dpath, fn)
        with open(prmt_fpath, 'rb') as fp:
            prmt = pickle.load(fp)
        T, _delta= [prmt.get(k) for k in ['T', '_delta']]
        new_row = [prefix, len(T), _delta]
        sol_fpath = opath.join(sol_dpath, 'sol_%s_%s.csv' % (prefix, aprc))
        log_fpath = opath.join(log_dpath, '%s_itr%s.csv' % (prefix, aprc))
        if opath.exists(sol_fpath):
            with open(sol_fpath) as r_csvfile:
                reader = csv.DictReader(r_csvfile)
                for row in reader:
                    objV, eliCpuTime = [row[cn] for cn in ['objV', 'eliCpuTime']]
            with open(opath.join(sol_dpath, 'sol_%s_%s.pkl' % (prefix, aprc)), 'rb') as fp:
                sol = pickle.load(fp)
            C, q_c = [sol.get(k) for k in ['C', 'q_c']]
            vehicles = [C[c] for c in range(len(C)) if q_c[c] > 0.5]
            handledTasks = list(chain(*vehicles))
            new_row += [eliCpuTime, objV,
                        len(vehicles), len(handledTasks),
                        len(handledTasks) / float(len(vehicles))]
        elif opath.exists(log_fpath):
            with open(log_fpath) as r_csvfile:
                reader = csv.DictReader(r_csvfile)
                for row in reader:
                    pass
                relObjV, eliCpuTime = [row[cn] for cn in ['relObjV', 'eliCpuTime']]
            new_row += [eliCpuTime, objV,
                        '-', '-',
                        '-']
        else:
            new_row += ['-', '-',
                        '-', '-',
                        '-']
        with open(sum_fpath, 'a') as w_csvfile:
            writer = csv.writer(w_csvfile, lineterminator='\n')
            writer.writerow(new_row)
    df = pd.read_csv(sum_fpath)
    df['seedNum'] = df.apply(lambda row: int(row['pn'].split('-')[-1][len('sn'):]), axis=1)
    df = df.sort_values(by=['seedNum'])
    df = df.drop(['seedNum'], axis=1)
    df.to_csv(sum_fpath, index=False)


def gen_problems(problem_dpath):
    if not opath.exists(problem_dpath):
        os.mkdir(problem_dpath)
        for dname in ['dplym', 'prmt']:
            os.mkdir(opath.join(problem_dpath, dname))
    summaryNV_dpath = opath.join(exp_dpath, '_summaryNV')
    dplym_dpath = reduce(opath.join, [summaryNV_dpath, 'dplym'])
    prmt_dpath = reduce(opath.join, [summaryNV_dpath, 'prmt'])
    PD_pairs = []
    for fn in os.listdir(prmt_dpath):
        if not fn.endswith('.pkl'): continue
        _, prefix = fn[:-len('.pkl')].split('_')
        with open(opath.join(prmt_dpath, fn), 'rb') as fp:
            prmt = pickle.load(fp)
        with open(opath.join(dplym_dpath, 'dplym_%s.pkl' % prefix), 'rb') as fp:
            flow_oridest, task_ppdp = pickle.load(fp)
        PD_pairs += task_ppdp
    #
    flowDplym, _, _mDP, _mTN, _dp, _fp, _ = prefix.split('-')
    NT = len(prmt['T'])
    dplym_dpath = reduce(opath.join, [problem_dpath, 'dplym'])
    prmt_dpath = reduce(opath.join, [problem_dpath, 'prmt'])

    for numTasks in np.arange(300, 301, 25):
        for seedNum in range(20):
            seed(seedNum)
            # numTasks = int(normalvariate(NT, NT * PER))
            task_ppdp = sample(PD_pairs, numTasks)
            problemName = '%s-nt%d-%s-%s-%s-%s-sn%d' % (flowDplym, numTasks, _mDP, _mTN, _dp, _fp, seedNum)
            with open(reduce(opath.join, [dplym_dpath, 'dplym_%s.pkl' % problemName]), 'wb') as fp:
                pickle.dump([flow_oridest, task_ppdp], fp)
            prmt = convert_prob2prmt(problemName, task_ppdp, numVeh=4)
            with open(opath.join(prmt_dpath, 'prmt_%s.pkl' % problemName), 'wb') as fp:
                pickle.dump(prmt, fp)


def summaryPA():
    summaryPA_dpath = opath.join(exp_dpath, '_summaryPA')
    rd_fpath = reduce(opath.join, [summaryPA_dpath, 'rawDataPA_VRP.csv'])
    with open(rd_fpath, 'w') as w_csvfile:
        writer = csv.writer(w_csvfile, lineterminator='\n')
        header = ['pn', 'numTasks', 'allowTT', 'cpuT',
                  'NHT', 'RHT', 'NV', 'XHNT',
                  'Revenue', 'Cost', 'OC', 'Profit']
        writer.writerow(header)
    prmt_dpath = reduce(opath.join, [summaryPA_dpath, 'prmt'])
    sol_dpath = reduce(opath.join, [summaryPA_dpath, 'sol'])
    aprc = 'CWL4'
    for fn in os.listdir(prmt_dpath):
        if not fn.endswith('.pkl'): continue
        _, prefix = fn[:-len('.pkl')].split('_')
        prmt_fpath = opath.join(prmt_dpath, fn)
        with open(prmt_fpath, 'rb') as fp:
            prmt = pickle.load(fp)
        T, _delta, cP, cC = [prmt.get(k) for k in ['T', '_delta', 'cP', 'cC']]
        new_row = [prefix, len(T), _delta]
        #
        sol_fpath = opath.join(sol_dpath, 'sol_%s_%s.csv' % (prefix, aprc))
        if not opath.exists(sol_fpath):
            new_row += ['-',
                        '-', '-', '-',
                        '-', '-', '-',
                        '-',
                        '-']
        else:
            with open(sol_fpath) as r_csvfile:
                reader = csv.DictReader(r_csvfile)
                for row in reader:
                    objV, eliCpuTime = [eval(row[cn]) for cn in ['objV', 'eliCpuTime']]
            with open(opath.join(sol_dpath, 'sol_%s_%s.pkl' % (prefix, aprc)), 'rb') as fp:
                sol = pickle.load(fp)
            C, q_c = [sol.get(k) for k in ['C', 'q_c']]
            vehicles = [C[c] for c in range(len(C)) if q_c[c] > 0.5]
            vn = len(vehicles)
            HNT = len(list(chain(*vehicles)))
            XHNT = len(T) - HNT
            assert abs(HNT * cP - objV) < 0.00001
            Revenue, Cost, OC = HNT * cP, len(vehicles) * cC, (len(T) - HNT) * cP
            Profit = Revenue - Cost - OC
            new_row += [eliCpuTime,
                        HNT, HNT / float(len(T)), vn, XHNT,
                        Revenue, Cost, OC,
                        Profit]
        with open(rd_fpath, 'a') as w_csvfile:
            writer = csv.writer(w_csvfile, lineterminator='\n')
            writer.writerow(new_row)
    df = pd.read_csv(rd_fpath)
    df['seedNum'] = df.apply(lambda row: int(row['pn'].split('-')[-1][len('sn'):]), axis=1)
    df = df.sort_values(by=['numTasks', 'seedNum'])
    df = df.drop(['seedNum'], axis=1)
    df.to_csv(rd_fpath, index=False)
    #
    sum_fpath = reduce(opath.join, [summaryPA_dpath, 'summaryPA_VRP.csv'])
    if not df[(df['cpuT'] == '-')].empty:
        df = df[(df['cpuT'] != '-')]
    df = df.drop(['pn'], axis=1)
    cols = list(df.columns[1:])
    df[cols] = df[cols].apply(pd.to_numeric)
    df = df.groupby(['numTasks']).mean().reset_index()
    df.to_csv(sum_fpath, index=False)


if __name__ == '__main__':
    gen_prmts(opath.join(exp_dpath, 'm100'), [5, 6, 7, 8])
    # run_experiments(2)
    # summaryNV()
    #
    # gen_problems(opath.join(exp_dpath, 'm102'))
    # summaryPA()
