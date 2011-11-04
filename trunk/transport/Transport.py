import numpy as np
from ase.transport.calculators import TransportCalculator

class TransportTB:
    """
    """
    def __init__(self,cal_LL,cal_SC,cal_LR,st_LL,st_SC,st_LR,vg_LL=0.0,vg_SC=0.0,vg_LR=0.0,txt=None):
        """
            This module is  used to  get the transport 
            properties the different unidimentional 
            system
        """
        self.st_LL=st_LL
        self.st_SC=st_SC
        self.st_LR=st_LR
        self.vg_LL=vg_LL
        self.vg_SC=vg_SC
        self.vg_LR=vg_LR
        self.cal_LL=cal_LL
        self.cal_SC=cal_SC
        self.cal_LR=cal_LR
        
    def get_transmmition(self,ene=np.arange(-1.0,1.0,0.01)):
        """
            This function return the Transmition function, 
            for a windows
        """
        self.cal_LL.solve_ground_state(self.st_LL)
        self.cal_SC.solve_ground_state(self.st_SC)
        self.cal_LR.solve_ground_state(self.st_LR)
        hLL,sLL,d,ds=self.cal_LL.ia.get_matrices()
        hSC,sSC,d,ds=self.cal_SC.ia.get_matrices()
        hLR,sLR,d,ds=self.cal_LR.ia.get_matrices()
        if self.cal_LL.get("SCC"):
            hL1=self.cal_LL.st.es.get_h1()
            hlp=hLL[0]+hL1*sLL[0]
        else:
            hlp=hLL[0]
        if self.cal_SC.get("SCC"):
            hSC1=self.cal_SC.st.es.get_h1()
            hsp=hSC[0]+hSC1*sSC[0]
        else:
            hsp=hSC[0]
        if self.cal_LR.get("SCC"):
            hsLR=self.cal_LR.st.es.get_h1()
            hlr=hLR[0]+hLR*sLR[0]
        else:
            h1r=hLR[0]
        self.tcalc = TransportCalculator(h=hsp,h1=hlp,h2=hlp,s=sSC[0],s1=sLL[0],s2=sLR[0],dos=False)
        self.tcalc.set(energies=ene)
        T_e = self.tcalc.get_transmission()
        return T_e
        
        