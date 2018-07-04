import scipy
import numpy,os
from sys import path
path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen')
import orbital_selection_fc as orb

from pyscf import gto,scf,cc,mp, symm, ao2mo, mp

def read_geometry(fname):
    inf = open(fname,'r')
    x   = inf.readlines()
    x   = x[2:]
    atm = []
    for z in x:
        y = z.split()
        rx,ry,rz = float(y[1]),float(y[2]),float(y[3])
        atm.append([y[0],numpy.asarray([rx,ry,rz])])
    return atm

#------------------------------------------
fname="/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_serial/test_riplinger_copy/DLPNO_CCSD_T-riplinger/diclophenac.xyz" #in files

""" diclophenac """
atoms = read_geometry(fname)
mol = gto.M(atom = atoms, basis='cc-pvdz') # read the geometry, in A
m = scf.RHF(mol)
m.kernel()

#mycc = cc.CCSD(m)
#mycc.kernel()

mymp2 = mp.MP2(m)
mymp2.kernel()

