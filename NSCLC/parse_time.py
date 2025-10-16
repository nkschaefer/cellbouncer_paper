#! /usr/bin/env python3
import sys
import os
import re

def main(args):
    
    fn = args[1]
    fnsplit = fn.split('_')
    
    prog = None
    experiment = None
    num = None
    
    if fnsplit[0] == "demux":
        if fnsplit[1] == 'vcf':
            prog = "CellBouncer"
        if fnsplit[2][0:2] == 'ds':
            experiment = "downsample"
            num = fnsplit[2].split('.time')[0][2:]
        elif fnsplit[2] == 'sc':
            experiment = 'ncell'
            num = fnsplit[3].split(".time")[0]
        elif fnsplit[2] == 'donor':
            experiment = 'nvar'
            num = fnsplit[-1].split(".time")[0]
    else:
        prog = fnsplit[0]
        if fnsplit[1] == "donor":
            experiment = 'nvar'
            num = fnsplit[-1].split('.time')[0]
        elif fnsplit[1][0:2] == 'ds':
            experiment = 'downsample'
            num = fnsplit[1].split('.time')[0][2:]
        elif fnsplit[1] == 'sc':
            experiment = 'ncell'
            num = fnsplit[2].split('.time')[0]
    
    f = open(args[1], 'r')
    time_sec = 0
    for line in f:
        line = line.rstrip()
        if line != "":
            dat = line.split('\t')
            if dat[0] == 'real':
                tsplit = dat[1].rstrip('s').split('m')
                tmin = int(tsplit[0])
                tsec = float(tsplit[1])
                time_sec = tmin * 60 + tsec
                break
    f.close()
    
    num = num.split('.time')[0]
    num = num.split('k')[0]
    print("{}\t{}\t{}\t{}".format(time_sec, prog, experiment, num))

if __name__ == '__main__':
    sys.exit(main(sys.argv))
