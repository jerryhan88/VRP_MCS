import os.path as opath
import os
import pickle
import shutil
#
from __path_organizer import exp_dpath
from problems import convert_prob2prmt
#
from _util_cython import gen_cFile


def gen_prmts(machine_dpath):
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
        prmt = convert_prob2prmt(problemName, task_ppdp)
        with open(opath.join(prmt_dpath, 'prmt_%s.pkl' % problemName), 'wb') as fp:
            pickle.dump(prmt, fp)


def run_experiments(machine_num):
    gen_cFile('CWL4')
    from CWL4 import run as CWL4_run
    machine_dpath = opath.join(exp_dpath, 'm%d' % machine_num)
    assert opath.exists(machine_dpath)
    gen_prmts(machine_dpath)
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


if __name__ == '__main__':
    # gen_prmts(opath.join(exp_dpath, 'm1'))
    run_experiments(2)