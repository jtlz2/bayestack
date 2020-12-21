'''
Eliab Malefahlo May 2019
This program computes the difference of log-evidence and draws them on a latex table
'''
import sys
import numpy
from utils import reportRelativeEvidences

chainz = []
red = [0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0]
'''
chainz.append(['101a','102a','202a_z_8x_2'])
chainz.append(['101b','102b','202b_z_8x_2'])
chainz.append(['101c','102c','202c_z_8x_2'])
chainz.append(['101d','102d','202d_z_8x_2'])
chainz.append(['101e','102e','202e_z_8x_2'])
chainz.append(['101f','102f','202f_z_8x_2'])
chainz.append(['101g','102g','202g_z_8x_2'])

chainz.append(['101h','102h','202h_z_8x_2'])
chainz.append(['101i','102i','202i_z_8x_2'])
chainz.append(['101j','102j','202j_z_8x_2'])
chain ='chains_200'
'''
chainz.append(['101a','102a','202a_z_8x_6_1'])
chainz.append(['101b','102b','202b_z_8x_6_1'])
chainz.append(['101c','102c','202c_z_8x_6_1'])
chainz.append(['101d','102d','202d_z_8x_6_1'])
chainz.append(['101e','102e','202e_z_8x_6_1'])
chainz.append(['101f','102f','202f_z_8x_6_1'])
chainz.append(['101g','102g','202g_z_8x_6_1'])

chainz.append(['101h','102h','202h_z_8x_6_1'])
chainz.append(['101i','102i','202i_z_8x_6_1'])
chainz.append(['101j','102j','202j_z_8x_6_1'])
chain ='chains_201'

#Epiing near bontvel station  (strong affordable shoes)

er =[]
z  = []

for i in range(len(chainz)):
    chains = [chain+a for a in chainz[i]]
    z_,er_ = reportRelativeEvidences(chains)
    z.append( z_)
    er.append(er_)

alpha=['A','B','C']
print z
print er
print len(z), len(er)

#for i in range(3):
#    print '%s & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$'%(alpha[i], z[0][i],er[0][i] , z[1][i],er[1][i] , z[2][i],er[2][i], z[3][i],er[3][i], z[4][i],er[4][i], z[5][i],er[5][i], z[6][i],er[6][i]) SDSS
print '\hline'
print 'Model & $%4.2f <z< %4.2f$ & $%4.2f <z< %4.2f$ & $%4.2f <z< %4.2f$ & $%4.2f <z< %4.2f$ & $%4.2f <z< %4.2f$  \hline'%(red[0],red[1], red[1],red[2],  red[2],red[3], red[3],red[4], red[4],red[5] )

for i in range(3):
    print '%s & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$'%(alpha[i], z[0][i],er[0][i] , z[1][i],er[1][i] , z[2][i],er[2][i], z[3][i],er[3][i], z[4][i],er[4][i])
#ev1 = myRelativeEvidences('%s %s %s'%(chains[0],chains[1],chains[2]))
print '\hline'
print ' & $%4.2f <z< %4.2f$ & $%4.2f <z< %4.2f$ & $%4.2f <z< %4.2f$ & $%4.2f <z< %4.2f$ & $%4.2f <z< %4.2f$'%(red[5],red[6], red[6],red[7], red[7],red[8],  red[8],red[9] ,red[9],red[10])
print '\hline'
for i in range(3):
    print '%s & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$ & ${%3.1f \pm %4.2f}$'%(alpha[i], z[5][i],er[5][i], z[6][i],er[6][i], z[7][i],er[7][i] , z[8][i],er[8][i], z[9][i],er[9][i])
    
    


