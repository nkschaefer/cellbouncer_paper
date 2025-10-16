#! /usr/bin/env python3
import sys
import os
import random
import gzip

def main(args):
    if len(args) < 4:
        print("ARGS: mtx features frac", file=sys.stderr)
        exit(1)

    mtx = args[1]
    features = args[2]
    frac = float(args[3])

    if frac <= 0 or frac >= 1:
        print("frac must be between 0 and 1", file=sys.stderr)
        exit(1)

    feature_inds = set([])
    fi = 1
    f = gzip.open(features, 'r')
    for line in f:
        line = line.decode().rstrip()
        dat = line.split('\t')
        if dat[2] == 'CRISPR Guide Capture':
            feature_inds.add(fi)
        fi += 1
    f.close()

    f = gzip.open(mtx, 'r')
    # Get tot
    first = True
    
    tot = 0
    newtot = 0
    
    fn_out = mtx.split('.')[0] + "_ds{}.matrix.mtx.gz".format(frac)
    f_out = gzip.open(fn_out, 'wt')
    
    hdr_out = []
    lines_out = []
    nfeature = 0
    nbc = 0

    for line in f:
        line = line.decode().rstrip()
        if line[0] == "%":
            hdr_out.append(line)
        else:
            if first:
                dat = line.split()
                nfeature = int(dat[0])
                nbc = int(dat[1])
                first = False
            else:
                dat = line.split()
                fi = int(dat[0])
                if fi not in feature_inds:
                    continue
                else:
                    count = float(dat[2])
                    tot += count
                    newcount = 0
                    for i in range(0, round(count)):
                        r = random.random()
                        if r <= frac:
                            newcount += 1
                    if newcount > 0:
                        lines_out.append('{} {} {}'.format(dat[0], dat[1], newcount))
                        newtot += newcount
    f.close()
    
    for h in hdr_out:
        print(h, file=f_out)
    print("{} {} {}".format(nfeature, nbc, len(lines_out)), file=f_out)
    for l in lines_out:
        print(l, file=f_out)

    f_out.close()

    print("Subsampled {} to {} total counts".format(tot, newtot), file=sys.stderr)



if __name__ == '__main__':
    sys.exit(main(sys.argv))
