from sys import path
path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen')
import orbital_selection_fc as orb
import numpy,os
from pyscf import gto,scf,cc,mp

""" penicillin"""

mol = gto.M(
    atom=[
    ['N',(3.17265,  1.15815, -0.09175)],
    ['C',(2.66167 , 0.72032 , 1.18601)],
    ['C',( 4.31931 , 0.59242 ,-0.73003)],
    ['C',(2.02252,  1.86922, -0.54680)],
    ['C',(1.37143,  1.52404,  0.79659)],
    ['S',(2.72625, -1.05563,  0.80065)],
    ['C',(4.01305, -0.91195, -0.52441)],
    ['C',(5.58297,  1.09423, -0.06535)],
    ['O',(1.80801,  2.36292, -1.62137)],
    ['N',(0.15715,  0.73759,  0.70095)],
    ['C',(5.25122, -1.72918, -0.12001)],
    ['C',(3.41769, -1.50152, -1.81857)],
    ['O',( 6.60623,  1.14077, -0.91855)],
    ['O',( 5.72538,  1.40990,  1.08931)],
    ['C',(-1.08932,  1.35001,  0.75816)],
    ['C',(-2.30230,  0.45820,  0.54941)],
    ['O',( -1.19855,  2.53493,  0.96288)],
    ['O',( -3.48875,  1.21403,  0.57063)],
    ['C',(-4.66939,  0.59150,  0.27339)],
    ['C',(-4.84065, -0.79240,  0.11956)],
    ['C',(-5.79523,  1.39165,  0.03916)],
    ['C',(-6.07568, -1.34753, -0.22401)],
    ['C',(-7.03670,  0.85454, -0.30482)],
    ['C',(-7.18253, -0.52580, -0.43612)],
    ['H',( 3.24354,  1.09074,  2.02120)],
    ['H',( 4.33865,  0.87909, -1.77554)],
    ['H',( 1.26605,  2.42501,  1.39138)],
    ['H',( 0.17381, -0.25857,  0.47675)],
    ['H',( 6.05024, -1.64196, -0.89101)],
    ['H',(  5.67754, -1.39089,  0.85176)],
    ['H',( 5.01118, -2.81229, -0.01401)],
    ['H',( 2.50304, -0.95210, -2.14173)],
    ['H',( 4.15186, -1.44541, -2.65467)],
    ['H',( 3.14138, -2.57427, -1.69700)],
    ['H',( 7.29069,  1.46408, -0.31004)],
    ['H',( -2.21049, -0.02915, -0.44909)],
    ['H',( -2.34192, -0.28647,  1.37775)],
    ['H',( -4.00164, -1.48999,  0.26950)],
    ['H',( -5.69703,  2.48656,  0.12872)],
    ['H',( -6.17811 ,-2.44045, -0.33185)],
    ['H',( -7.89945 , 1.51981, -0.47737)],
    ['H',( -8.15811, -0.96111, -0.71027)]],
    basis ='cc-pvdz')

m = scf.RHF(mol)
m.kernel()

#mycc = cc.CCSD(m)
#mycc.kernel()

mymp2 = mp.MP2(m)
mymp2.kernel()



