import numpy as np
import csv

def getColumn(filename, column, delimiter=',', skipinitialspace=False, skipheader=True, magnaprobe=False):
    results = csv.reader(open(filename),delimiter=delimiter,skipinitialspace=skipinitialspace)
    if skipheader==True:
        next(results, None)
    if magnaprobe==True: #has 4 headers
        next(results, None)
        next(results, None)
        next(results, None)
    return [result[column] for result in results]
