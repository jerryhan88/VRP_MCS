

def run(prmt, pd_inputs):
    t_ij = prmt['t_ij']
    seq0, i0 = list(map(pd_inputs.get, ['seq0', 'i0']))
    #
    iP, iM = 'p%d' % i0, 'd%d' % i0
    travelTime, seq1 = 1e400, None
    assert len(seq0) >= 2
    for i in range(1, len(seq0) - 1):
        for j in range(i, len(seq0) - 1):
            seq = seq0[:]
            seq.insert(i, iP)
            seq.insert(j + 1, iM)
            detourTime = calc_travelTime(seq, t_ij)
            if detourTime < travelTime:
                travelTime, seq1 = detourTime, seq
    return travelTime, seq1


def calc_travelTime(seq, t_ij):
    travelTime = 0.0
    for i in range(len(seq) - 1):
        travelTime += t_ij[seq[i], seq[i + 1]]
    return travelTime


if __name__ == '__main__':
    pass
