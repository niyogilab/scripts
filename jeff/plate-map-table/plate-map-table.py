#!/usr/bin/env python

'''@usageTxt@'''

from csv import reader, writer
from sys import stdout

def readcsv(filename):
    with open(filename, 'rb') as csv:
        return [[c for c in r] for r in reader(csv)]

def tidycsv(mapcsv, colname):
    rows = [row[0]  for row in mapcsv[1:]]
    vals = [row[1:] for row in mapcsv[1:]]
    cols = mapcsv[0][1:]
    cols2 = ['well', colname]
    vals2 = []
    for r in range(len(rows)):
        for c in range(len(cols)):
            vals2.append([rows[r]+cols[c], vals[r][c]])
    return [cols2] + vals2

def writecsv(rows, filename):
    with open(filename, 'wb') as csv:
        writer(csv).writerows(rows)

def main(args):
    csv = readcsv(args['<incsv>'])
    csv = tidycsv(csv, args['<colname>'])
    writecsv(csv, args['<outcsv>'])

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__, version='@name@')
    main(args)
