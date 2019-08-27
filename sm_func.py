import numpy as np
import csv

def getColumn(filename, column, delimiter=',', skipinitialspace=False, skipheader=True):
    results = csv.reader(open(filename),delimiter=delimiter,skipinitialspace=skipinitialspace)
    if skipheader==True:
        next(results, None)
    return [result[column] for result in results]
