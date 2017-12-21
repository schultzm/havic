
from rpy2.rinterface import RRuntimeError
from rpy2.robjects.packages import importr
utils = importr('utils')
contriburl = 'https://cran.ms.unimelb.edu.aux/'

#move this to class
def importr_tryhard(packname):
    try:
        rpack = importr(packname)
    except RRuntimeError:
        utils.install_packages(packname, contriburl = contriburl)
        rpack = importr(packname)
    return rpack