#####################################################################
#   Program  for Calculed  electronic Structur
#              J. D. Correa
#                 2011
#####################################################################


from ase import *
from hotbit import *
from box.interpolation import interpolate_path
import numpy  as np
from Transport import TransportTB
#--------------------------------------------------------------------
# setup  input strures: lead_left, conductor, lead_right
#--------------------------------------------------------------------
nt_con = read('nt100_lead.xyz',index=-1,format='xyz') 		#read atomic position  of the conductor
nt_leadL = read('nt100_lead.xyz',index=-1,format='xyz')		#read atomic position  of the left lead
nt_leadR = read('nt100_lead.xyz',index=-1,format='xyz')		#read atomic position  of the right lead
#--------------------------------------------------------------------
# control var to  include vgap  and perpendicular electric field
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#   Configure hotbit calculation options
#--------------------------------------------------------------------
tab={'CC':'H_H.par'}
elem={'C':'H.elm'}
#,tables=tab,elements=elem
cal = Hotbit(SCC=True, charge_density= 'Gaussian',width=0.01, rs='physical_k',maxiter=50,txt='LL.out',kpts=(1,1,1),gamma_cut=7.0)
cal2 = Hotbit(SCC=True, charge_density= 'Gaussian',width=0.01, rs='physical_k',maxiter=50,txt='SS.out',kpts=(1,1,1),gamma_cut=7.0)
cal3 = Hotbit(SCC=True, charge_density= 'Gaussian',width=0.01, rs='physical_k',maxiter=50,txt='LR.out',kpts=(1,1,1),gamma_cut=7.0)
tb=TransportTB(cal_LL=cal2,cal_SC=cal2,cal_LR=cal3,st_LL=nt_leadL,st_SC=nt_con,st_LR=nt_leadR)
e=np.arange(-0.4,0.4,0.001)
T_e=tb.get_transmmition(ene=e)
f=open("test_transmission.dat", "w")
for i  in range(np.size(T_e)):
    f.write(str(e[i]*27.211)+"\t"+str(T_e[i])+"\n")
