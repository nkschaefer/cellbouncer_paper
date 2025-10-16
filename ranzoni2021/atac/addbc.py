#! /usr/bin/env python3
import os
import sys
import gzip

def proc_file(fn1, fn1out, bc):
    f = gzip.open(fn1, 'r')
    fout = gzip.open(fn1out, 'wb')

    for line in f:
        line = line.decode().rstrip()
        if line[0] == '@':
            line = line.split(' ')[0]
            line += ' CB:Z:{}'.format(bc)
        line += "\n"
        fout.write(line.encode())
         
    f.close()
    fout.close()

def main(args):

    fn1 = args[1]
    fn2 = args[2]
    outd = args[3]
    bc = args[4]
    
    if outd[-1] != '/':
        outd += '/'
    fn1out = outd + fn1.split('/')[-1]
    fn2out = outd + fn2.split('/')[-1]
    
    proc_file(fn1, fn1out, bc)
    proc_file(fn2, fn2out, bc)

if __name__ == '__main__':
    sys.exit(main(sys.argv))

