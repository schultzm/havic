from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
utils = importr('utils')

def importr_tryhard(packname):
    contriburl = 'http://cran.stat.ucla.edu/'
    try:
        rpack = importr(packname, on_conflict="warn")
    except RRuntimeError:
        rpack = []
    return rpack
