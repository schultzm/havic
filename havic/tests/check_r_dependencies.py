from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeWarning
utils = importr('utils')

def importr_tryhard(packname):
    try:
        rpack = importr(packname, on_conflict="warn")
    except RRuntimeWarning:
        rpack = []
    return rpack
