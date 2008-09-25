from repulsion import RepulsivePotential
from box import mix
from box.interpolation import Function
from hotbit.fortran.slako import fast_slako_transformations
import numpy as nu
from hotbit import auxil
from box.mix import kronecker
from box.interpolation import MultipleSplineFunction
dot=nu.dot
array=nu.array

aux={'s':0,'p':1,'d':2}
itypes={'ss':['s'],'sp':['s'],'ps':['s'],'sd':['s'],'ds':['s'],
        'pp':['s','p'],'pd':['s','p'],'dp':['s','p'],'dd':['s','p','d']}
integrals={'dds':0,'ddp':1,'ddd':2,'pds':3,'pdp':4,'pps':5,'ppp':6,\
           'sds':7,'sps':8,'sss':9,'dps':10,'dpp':11,'dss':12,'pss':13}                






class Interactions:
    
    def __init__(self,calc,timer,elements,tables):
        """
        Set up the input files for interactions.
        
        If tables==None, use default tables.
        If tables['others']=='default' use default interactions for the ones not given in tables.
        
        """ 
        from os import environ 
        from os.path import isfile
        present=elements.get_present()
        default_dir=environ.get('HOTBIT_PARAMETERS')
        files={}
        if tables!=None:
            tables=tables.copy()
        
        if tables==None or ('others' in tables and tables['others']=='default'):
            if tables!=None and 'others' in tables:
                tables.pop('others')
            for e1 in present:
                for e2 in present:
                    # select file which is present (e.g. H_C.par or C_H.par)
                    file12='./%s_%s.par' %(e1,e2)
                    file21='./%s_%s.par' %(e2,e1)
                    default12='%s/%s_%s.par' %(default_dir,e1,e2)
                    default21='%s/%s_%s.par' %(default_dir,e2,e1)
                    if isfile(file12):
                        file=file12
                    elif isfile(file21):
                        file=file21
                    elif isfile(default12):
                        file=default12
                    elif isfile(default21):
                        file=default21
                    else:         
                        file='not found in default place'
                    files[e1+e2]=files[e2+e1]=file            
        
        if tables!=None:
            for table in tables:
                if table=='others':
                    continue
                # search table first from present dir, then from default
                file=tables[table]
                default='%s/%s' %(default_dir,file)
                e1,e2=auxil.separate_symbols(table)
                if isfile(file):   
                    pass
                elif isfile(default):
                    file=default
                else:
                    file='not found in custom place'
                files[e1+e2]=files[e2+e1]=file
                
        for file in files:
            if not isfile(files[file]):
                raise AssertionError('Interaction for %s %s.' %(file,files[file]))
        self.files=files
        self.calc=calc
        self.present=present
        self.timer=timer
        self.el=elements
        self.max_cut=0.0 # maximum interaction range in Bohrs
        self.read_tables()
        self.check_box_size()
        self.first=True
    
    def get_files(self):
        """ Return the list of Slater-Koster table files."""
        return self.files  
         
    def get_cutoffs(self):
        """ Return element pair cutoff dictionary. """
        return self.cut         
    
    def check_box_size(self):
        """ Check that box is large enough. """
        gamma_cut, SCC=self.calc.get('gamma_cut') ,self.calc.get('SCC')
        pbc=self.el.atoms.get_pbc()
        for periodic,length in zip(self.el.get_box_lengths(),pbc):
            for ia in self.cut:
                if self.cut[ia]>length/2 and not periodic:
                    raise AssertionError('Too small box (one size %.2f, largest cut %.2f)' %(length,self.cut[ia]))
                if SCC and periodic and (gamma_cut!=None and gamma_cut<length/4):
                    raise AssertionError('gamma_cut should be small enough compared to %s Bohr' %length)
                     
    def greetings(self):
        """ Return the interaction specifications """
        txt='Interactions:\n'
        shown=[]
        for ia in self.files:
            file=self.files[ia]
            if file in shown: continue
            txt+='  %s in %s\n' %(ia,file)
            doc=mix.find_value(file,'slako_comment',fmt='strings',default=['no slako doc'])
            for line in doc:
                txt+='    *'+line.lstrip()+'\n'     
            shown.append(file)                
        return txt
             
    def read_tables(self):
        """ 
        Read par-file tables. Files are tabulated with |ket> having >= angular momentum,
        so one has to copy values for interaction the other way round.
        """
        self.h={}
        self.s={}
        self.cut={}
        self.maxh={}          
                
        for si in self.present:
            for sj in self.present:
                try:
                    table_ij=mix.find_value(self.files[si+sj],'%s_%s_table' %(si,sj),fmt='matrix')
                    table_ji=mix.find_value(self.files[sj+si],'%s_%s_table' %(sj,si),fmt='matrix')
                except KeyError:
                    raise KeyError('Interaction file for %s-%s or %s-%s not found.' %(si,sj,sj,si))
                self.cut[si+sj]=table_ij[-1,0]
                self.max_cut=max(self.max_cut,self.cut[si+sj])
                self.maxh[si+sj]=max( [max(nu.abs(table_ij[:,i])) for i in range(1,11)] )
                
                ei, ej=self.el.elements[si], self.el.elements[sj]
                valence_i, valence_j=ei.get_valence_orbitals(), ej.get_valence_orbitals()
                
                pair=si+sj
                self.h[pair]=MultipleSplineFunction(table_ij[:,0])
                self.s[pair]=MultipleSplineFunction(table_ij[:,0])
                for vi in valence_i:
                    for vj in valence_j:
                        li, lj=vi[1], vj[1]
                        # for given valence orbitals, go through all possible integral types (sigma,pi,...)
                        for itype in itypes[li+lj]:
                            table='%s(%s)-%s(%s)-%s' %(si,li,sj,lj,itype)
                            short='%s%s%s' %(li,lj,itype)
                            # self.h['C(p)-H(s)-sigma']=...
                            if short[0:2]=='ps' or short[0:2]=='ds' or short[0:2]=='dp':
                                # this is tabulated in other table; switch order -> parity factor
                                parity=(-1)**( aux[li]+aux[lj] )
                                index=integrals[short[1]+short[0]+short[2]]
                                self.h[pair].add_function(table_ji[:,index+1]*parity,table,integrals[short])
                                self.s[pair].add_function(table_ji[:,index+11]*parity,table,integrals[short])
                            else:
                                index=integrals[short]                                
                                self.h[pair].add_function(table_ij[:,index+1],table,integrals[short])
                                self.s[pair].add_function(table_ij[:,index+11],table,integrals[short])
                                
        
    def plot_table(self,e1,e2,der=0):
        """ Plot SlaKo table for given elements. """
        import pylab as pl
        R=nu.linspace(0,self.cut[e1+e2],100)
        raise NotImplementedError('...h[pair]')
        for i in range(10):
            pl.subplot(5,2,i+1)
            pl.plot(R,[self.h[e1+e2][i](r,der=der) for r in R],label='H')
            pl.plot(R,[self.s[e1+e2][i](r,der=der) for r in R],label='S')
            pl.xlim(0,R[-1])
            pl.ylim(-self.maxh[e1+e2],self.maxh[e1+e2])
            pl.title('%s-%s integral %s (der=%i)' %(e1,e2,self.integrals[i],der))
            pl.xlabel('r (Bohr)')
            pl.ylabel('H (Hartree) and S')
        pl.show()
        
    def construct_matrices(self):
        """ Hamiltonian and overlap matrices. """
        start=self.timer.start
        stop=self.timer.stop
        start('matrix construction')
        orbs=self.el.orbitals()
        norb=len(orbs)
        self.H0=nu.zeros((norb,norb))
        self.dH0, self.dS=nu.zeros((norb,norb,3)),nu.zeros((norb,norb,3))
        self.S=nu.identity(norb)            
        
        orbitals=[[orb['orbital'] for orb in self.el.orbitals(i)] for i in range(len(self.el))]
        orbindex=[self.el.orbitals(i,indices=True) for i in range(len(self.el))]
        
        self.el.set_cutoffs(self.cut)
        nonzero=0
        for i,j,si,sj,dist,r,rhat in self.el.get_ia_atom_pairs(['i','j','si','sj','dist','r','rhat']):
            if i==j:
                for orb in self.el.orbitals(i):
                    self.H0[orb['index'],orb['index']]=orb['energy']
                    nonzero+=1
            else:                    
                # make interpolated mels, also derivatives
                h, s, dh, ds=nu.zeros((14,)), nu.zeros((14,)), nu.zeros((14,3)), nu.zeros((14,3))               
                pair=si+sj
                               
                start('mel splint')                                                                
                hi, dhi=self.h[pair](dist)
                si, dsi=self.s[pair](dist)
                stop('mel splint')                                                                    
                
                start('setup pre-h')         
                indices=self.h[pair].get_indices()                                                                      
                h[indices], s[indices]=hi, si                
                dh[indices], ds[indices]=nu.outer(dhi,rhat), nu.outer(dsi,rhat)
                stop('setup pre-h')                                       
     
                start('fortran slako')
                obsi, obsj=orbindex[i], orbindex[j]
                noi, noj=len(obsi), len(obsj)
                #ht,st,dht,dst=slako_transformations(rhat,dist,noi,noj,h,s,dh,ds)
                ht, st, dht, dst=fast_slako_transformations(rhat,dist,noi,noj,h,s,dh,ds)
                
                i1, i2, j1, j2=obsi[0], obsi[-1]+1, obsj[0], obsj[-1]+1
                self.H0[i1:i2,j1:j2]=ht
                self.S[i1:i2,j1:j2]=st
                self.dH0[i1:i2,j1:j2,:]=dht
                self.dS[i1:i2,j1:j2,:]=dst 
                
                # symmetrize (antisymmetrize) H and S (dH and dS)
                self.H0[j1:j2,i1:i2]=ht.transpose()
                self.S[j1:j2,i1:i2]=st.transpose()
                self.dH0[j1:j2,i1:i2,:]=-dht.transpose((1,0,2))
                self.dS[j1:j2,i1:i2,:]=-dst.transpose((1,0,2))                 
                
                if self.first:
                    nonzero+=sum( abs(ht.flatten())>1E-20 )*2
                stop('fortran slako')
                                    
        if self.first:
            self.calc.out('Hamiltonian matrix is %.3f %% filled on first calculation.' %(nonzero*100.0/norb**2) )                     
            self.first=False                    
                               
        stop('matrix construction')
            
    def get_matrices(self):
        self.construct_matrices()
        return self.H0,self.S,self.dH0,self.dS
    
    def get_cut(self):
        """ Maximum cutoff. """
        return self.max_cut
    
def simple_table_notation(table):
    a,b,i=table[2:].split('-')
    return a[-2]+b[-2]+i[0]                
            
            
s3=nu.sqrt(3.0) 

def slako_transformations(rhat,dist,noi,noj,h,s,dh,ds):
    """ 
    Apply Slater-Koster transformation rules to orbitals iorbs and orbitals jorbs,
    where rhat is vector i->j and table gives the values for given tabulated
    matrix elements. Convention: orbital name starts with s,p,d,...
    """    
    l,m,n=rhat          
    ll,mm,nn=rhat**2 
    dl=(nu.array([1,0,0])-l*rhat)/dist
    dm=(nu.array([0,1,0])-m*rhat)/dist
    dn=(nu.array([0,0,1])-n*rhat)/dist
    dll, dmm, dnn = 2*l*dl, 2*m*dm, 2*n*dn 
    
    mat=nu.zeros((9,9,14))
    ind=nu.zeros((9,9,14),dtype=int)
    der=nu.zeros((9,9,14,3))
    cnt=nu.zeros((9,9),dtype=int)+1
    cnt[1:,1:]=2
    cnt[4:,4:]=3
    
    ht=nu.zeros((noi,noj))
    st=nu.zeros((noi,noj))
    dht=nu.zeros((noi,noj,3))
    dst=nu.zeros((noi,noj,3))
    mxorb=max(noi,noj)
    
    mat[0,0,0]=1  #ss
    der[0,0,0]=0
    ind[0,0,0]=9
    
    if mxorb>=2:   #sp
        mat[0,1,0]=l
        der[0,1,0,:]=dl
        ind[0,1,0]=8
            
        mat[0,2,0]=m
        der[0,2,0,:]=dm                    
        ind[0,2,0]=8
                
        mat[0,3,0]=n
        der[0,3,0,:]=dn 
        ind[0,3,0]=8
    
    if mxorb>=5:   #sd            
        mat[0,4,0]=s3*l*m
        der[0,4,0,:]=s3*(dl*m+l*dm)
        ind[0,4,0]=7
                
        mat[0,5,0]=s3*m*n
        der[0,5,0,:]=s3*(dm*n+m*dn)                
        ind[0,5,0]=7
                
        mat[0,6,0]=s3*n*l
        der[0,6,0,:]=s3*(dn*l+n*dl)            
        ind[0,6,0]=7
                
        mat[0,7,0]=0.5*s3*(ll-mm)
        der[0,7,0,:]=0.5*s3*(dll-dmm)            
        ind[0,7,0]=7
                
        mat[0,8,0]=nn-0.5*(ll+mm)
        der[0,8,0,:]=dnn-0.5*(dll+dmm)
        ind[0,8,0]=7
               
    if mxorb>=2: #pp
        mat[1,1,0:2]=[ll, 1-ll]
        der[1,1,0:2,:]=[dll, -dll]
        ind[1,1,0:2]=[5,6]
            
        mat[1,2,0:2]=[l*m, -l*m]
        der[1,2,0:2,:]=[dl*m+l*dm, -(dl*m+l*dm)]
        ind[1,2,0:2]=[5,6]
            
        mat[1,3,0:2]=[l*n, -l*n]
        der[1,3,0:2,:]=[dl*n+l*dn, -(dl*n+l*dn)]
        ind[1,3,0:2]=[5,6]
        
    if mxorb>=5: #pd                
        mat[1,4,0:2]=[s3*ll*m, m*(1-2*ll)]
        der[1,4,0:2,:]=[s3*(dll*m+ll*dm), dm*(1-2*ll)+m*(-2*dll)]
        ind[1,4,0:2]=[3,4]
                
        mat[1,5,0:2]=[s3*l*m*n, -2*l*m*n]
        der[1,5,0:2,:]=[s3*(dl*m*n+l*dm*n+l*m*dn), -2*(dl*m*n+l*dm*n+l*m*dn)]
        ind[1,5,0:2]=[3,4]
                
        mat[1,6,0:2]=[s3*ll*n, n*(1-2*ll)]
        der[1,6,0:2,:]=[s3*(dll*n+ll*dn), dn*(1-2*ll)+n*(-2*dll)]
        ind[1,6,0:2]=[3,4]
                
        mat[1,7,0:2]=[0.5*s3*l*(ll-mm), l*(1-ll+mm)]
        der[1,7,0:2,:]=[0.5*s3*(dl*(ll-mm)+l*(dll-dmm)), dl*(1-ll+mm)+l*(-dll+dmm)]
        ind[1,7,0:2]=[3,4]
                
        mat[1,8,0:2]=[l*(nn-0.5*(ll+mm)), -s3*l*nn]
        der[1,8,0:2,:]=[dl*(nn-0.5*(ll+mm))+l*(dnn-0.5*(dll+dmm)), -s3*(dl*nn+l*dnn)]
        ind[1,8,0:2]=[3,4]
            
    if mxorb>=2:            
        mat[2,2,0:2]=[mm, 1-mm]
        der[2,2,0:2,:]=[dmm, -dmm]
        ind[2,2,0:2]=[5,6]
                
        mat[2,3,0:2]=[m*n, -m*n]
        der[2,3,0:2,:]=[dm*n+m*dn, -(dm*n+m*dn)]
        ind[2,3,0:2]=[5,6]
            
    if mxorb>=5:            
        mat[2,4,0:2]=[s3*mm*l, l*(1-2*mm)]
        der[2,4,0:2,:]=[s3*(dmm*l+mm*dl), dl*(1-2*mm)+l*(-2*dmm)]
        ind[2,4,0:2]=[3,4]
                
        mat[2,5,0:2]=[s3*mm*n, n*(1-2*mm)]
        der[2,5,0:2,:]=[s3*(dmm*n+mm*dn), dn*(1-2*mm)+n*(-2*dmm)]
        ind[2,5,0:2]=[3,4]
                
        mat[2,6,0:2]=[s3*m*n*l, -2*m*n*l]
        der[2,6,0:2,:]=[s3*(dm*n*l+m*dn*l+m*n*dl), -2*(dm*n*l+m*dn*l+m*n*dl)]
        ind[2,6,0:2]=[3,4]
                
        mat[2,7,0:2]=[0.5*s3*m*(ll-mm), -m*(1+ll-mm)]
        der[2,7,0:2,:]=[0.5*s3*(dm*(ll-mm)+m*(dll-dmm)), -(dm*(1+ll-mm)+m*(dll-dmm))]
        ind[2,7,0:2]=[3,4]
                
        mat[2,8,0:2]=[m*(nn-0.5*(ll+mm)), -s3*m*nn]
        der[2,8,0:2,:]=[dm*(nn-0.5*(ll+mm))+m*(dnn-0.5*(dll+dmm)), -s3*(dm*nn+m*dnn)]
        ind[2,8,0:2]=[3,4]
            
    if mxorb>=2:
        mat[3,3,0:2]=[nn, 1-nn]
        der[3,3,0:2,:]=[dnn, -dnn]
        ind[3,3,0:2]=[5,6]
            
    if mxorb>=5:
        mat[3,4,0:2]=[s3*l*m*n, -2*m*n*l]
        der[3,4,0:2,:]=[s3*(dl*m*n+l*dm*n+l*m*dn), -2*(dm*n*l+m*dn*l+m*n*dl)]
        ind[3,4,0:2]=[3,4]
                
        mat[3,5,0:2]=[s3*nn*m, m*(1-2*nn)]
        der[3,5,0:2,:]=[s3*(dnn*m+nn*dm), dm*(1-2*nn)+m*(-2*dnn)]
        ind[3,5,0:2]=[3,4]
                
        mat[3,6,0:2]=[s3*nn*l, l*(1-2*nn)]
        der[3,6,0:2,:]=[s3*(dnn*l+nn*dl), dl*(1-2*nn)+l*(-2*dnn)]
        ind[3,6,0:2]=[3,4]
                
        mat[3,7,0:2]=[0.5*s3*n*(ll-mm), -n*(ll-mm)]
        der[3,7,0:2,:]=[0.5*s3*(dn*(ll-mm)+n*(dll-dmm)), -(dn*(ll-mm)+n*(dll-dmm))]
        ind[3,7,0:2]=[3,4]
                
        mat[3,8,0:2]=[n*(nn-0.5*(ll+mm)), s3*n*(ll+mm)]
        der[3,8,0:2,:]=[dn*(nn-0.5*(ll+mm))+n*(dnn-0.5*(dll+dmm)), s3*(dn*(ll+mm)+n*(dll+dmm))]
        ind[3,8,0:2]=[3,4]
    
    if mxorb>=5:
        mat[4,4,0:3]=[3*ll*mm, ll+mm-4*ll*mm, nn+ll*mm]
        der[4,4,0:3,:]=[3*(dll*mm+ll*dmm), dll+dmm-4*(dll*mm+ll*dmm), dnn+(dll*mm+ll*dmm)]
        ind[4,4,0:3]=[0,1,2]
            
        mat[4,5,0:3]= [3*l*mm*n, l*n*(1-4*mm), l*n*(mm-1)]
        der[4,5,0:3,:]= [3*(dl*mm*n+l*dmm*n+l*mm*dn), dl*n*(1-4*mm)+l*dn*(1-4*mm)+l*n*(-4*dmm), dl*n*(mm-1)+l*dn*(mm-1)+l*n*(dmm)]
        ind[4,5,0:3]=[0,1,2]
            
        mat[4,6,0:3]=[3*ll*m*n, m*n*(1-4*ll), m*n*(ll-1)]
        der[4,6,0:3,:]=[3*(dll*m*n+ll*dm*n+ll*m*dn), dm*n*(1-4*ll)+m*dn*(1-4*ll)+m*n*(-4*dll), dm*n*(ll-1)+m*dn*(ll-1)+m*n*(dll)]
        ind[4,6,0:3]=[0,1,2]
            
        mat[4,7,0:3]=[1.5*l*m*(ll-mm), 2*l*m*(mm-ll), 0.5*l*m*(ll-mm)]
        der[4,7,0:3,:]=[1.5*(dl*m*(ll-mm)+l*dm*(ll-mm)+l*m*(dll-dmm)),\
                    2*(dl*m*(mm-ll)+l*dm*(mm-ll)+l*m*(dmm-dll)),\
                    0.5*(dl*m*(ll-mm)+l*dm*(ll-mm)+l*m*(dll-dmm))]
        ind[4,7,0:3]=[0,1,2]                    
            
        mat[4,8,0:3]=[s3*l*m*(nn-0.5*(ll+mm)), - 2*s3*l*m*nn, 0.5*s3*l*m*(1+nn)]
        der[4,8,0:3,:]=[s3*( dl*m*(nn-0.5*(ll+mm))+l*dm*(nn-0.5*(ll+mm))+l*m*(dnn-0.5*(dll+dmm)) ),\
                    -2*s3*(dl*m*nn+l*dm*nn+l*m*dnn),\
                    0.5*s3*( dl*m*(1+nn)+l*dm*(1+nn)+l*m*(dnn) )]
        ind[4,8,0:3]=[0,1,2]                    
            
        mat[5,5,0:3]=[3*mm*nn,  (mm+nn-4*mm*nn), (ll+mm*nn)]
        der[5,5,0:3,:]=[3*(dmm*nn+mm*dnn), (dmm+dnn-4*(dmm*nn+mm*dnn)),  (dll+dmm*nn+mm*dnn)]            
        ind[5,5,0:3]=[0,1,2]
            
        mat[5,6,0:3]=[3*m*nn*l, m*l*(1-4*nn), m*l*(nn-1)]
        der[5,6,0:3,:]=[3*(dm*nn*l+m*dnn*l+m*nn*dl),\
                    dm*l*(1-4*nn)+m*dl*(1-4*nn)+m*l*(-4*dnn),\
                    dm*l*(nn-1)+m*dl*(nn-1)+m*l*(dnn)]              
        ind[5,6,0:3]=[0,1,2]                    
            
        mat[5,7,0:3]=[1.5*m*n*(ll-mm), - m*n*(1+2*(ll-mm)), m*n*(1+0.5*(ll-mm))]
        der[5,7,0:3,:]=[1.5*( dm*n*(ll-mm)+m*dn*(ll-mm)+m*n*(dll-dmm) ),\
                    - ( dm*n*(1+2*(ll-mm))+m*dn*(1+2*(ll-mm))+m*n*(2*dll-2*dmm) ),\
                    dm*n*(1+0.5*(ll-mm))+m*dn*(1+0.5*(ll-mm))+m*n*(0.5*(dll-dmm))]            
        ind[5,7,0:3]=[0,1,2]                    
            
        mat[5,8,0:3]=[s3*m*n*(nn-0.5*(ll+mm)), s3*m*n*(ll+mm-nn), -0.5*s3*m*n*(ll+mm)]
        der[5,8,0:3,:]=[s3*( dm*n*(nn-0.5*(ll+mm)) + m*dn*(nn-0.5*(ll+mm))+m*n*(dnn-0.5*(dll+dmm)) ),\
                    s3*( dm*n*(ll+mm-nn)+m*dn*(ll+mm-nn)+m*n*(dll+dmm-dnn) ),\
                    - 0.5*s3*( dm*n*(ll+mm)+m*dn*(ll+mm)+m*n*(dll+dmm) )]                   
        ind[5,8,0:3]=[0,1,2]                    
            
        mat[6,6,0:3]=[3*nn*ll, (nn+ll-4*nn*ll), (mm+nn*ll)]
        der[6,6,0:3,:]=[3*(dnn*ll+nn*dll), dnn+dll-4*(dnn*ll+nn*dll), (dmm+dnn*ll+nn*dll)]            
        ind[6,6,0:3]=[0,1,2]
            
        mat[6,7,0:3]=[1.5*n*l*(ll-mm), n*l*(1-2*(ll-mm)), - n*l*(1-0.5*(ll-mm))]
        der[6,7,0:3,:]=[1.5*( dn*l*(ll-mm)+n*dl*(ll-mm)+n*l*(dll-dmm) ),\
                    dn*l*(1-2*(ll-mm))+n*dl*(1-2*(ll-mm))+n*l*(-2*(dll-dmm)),\
                    -( dn*l*(1-0.5*(ll-mm))+n*dl*(1-0.5*(ll-mm))+n*l*(-0.5*(dll-dmm)) )]
        ind[6,7,0:3]=[0,1,2]                    
                
        mat[6,8,0:3]=[s3*l*n*(nn-0.5*(ll+mm)), s3*l*n*(ll+mm-nn), - 0.5*s3*l*n*(ll+mm)]
        der[6,8,0:3,:]=[s3*( dl*n*(nn-0.5*(ll+mm))+l*dn*(nn-0.5*(ll+mm))+l*n*(dnn-0.5*(dll+dmm)) ),\
                    s3*( dl*n*(ll+mm-nn)+l*dn*(ll+mm-nn)+l*n*(dll+dmm-dnn) ),\
                    - 0.5*s3*( dl*n*(ll+mm)+l*dn*(ll+mm)+l*n*(dll+dmm) )]             
        ind[6,8,0:3]=[0,1,2]                    
    
        mat[7,7,0:3]=[0.75*(ll-mm)**2, (ll+mm-(ll-mm)**2), (nn+0.25*(ll-mm)**2)]
        der[7,7,0:3,:]=[0.75*2*(ll-mm)*(dll-dmm), (dll+dmm-2*(ll-mm)*(dll-dmm)), (dnn+0.25*2*(ll-mm)*(dll-dmm))]                
        ind[7,7,0:3]=[0,1,2]        
    
        mat[7,8,0:3]=[0.5*s3*(ll-mm)*(nn-0.5*(ll+mm)), s3*nn*(mm-ll), 0.25*s3*(1+nn)*(ll-mm)]
        der[7,8,0:3,:]=[0.5*s3*( (dll-dmm)*(nn-0.5*(ll+mm))+(ll-mm)*(dnn-0.5*(dll+dmm)) ),\
                    s3*( dnn*(mm-ll)+nn*(dmm-dll) ),\
                    0.25*s3*( dnn*(ll-mm)+(1+nn)*(dll-dmm) )]            
        ind[7,8,0:3]=[0,1,2]
                            
        mat[8,8,0:3]=[(nn-0.5*(ll+mm))**2, 3*nn*(ll+mm), 0.75*(ll+mm)**2]
        der[8,8,0:3,:]=[2*(nn-0.5*(ll+mm))*(dnn-0.5*(dll+dmm)),\
                    3*( dnn*(ll+mm)+nn*(dll+dmm) ),\
                    0.75*2*(ll+mm)*(dll+dmm)]   
        ind[8,8,0:3]=[0,1,2]                    
                
    # use the same rules for orbitals when they are reversed (pd ->dp)...
    for a in range(9):
        for b in range(a+1,9):
            mat[b,a,:]=mat[a,b,:]
            der[b,a,:,:]=der[a,b,:,:]
            ind[b,a,:]=ind[a,b,:]
            
    # ...but use different indices from table            
    #pd 3:5-->10:12 
    #sd 7->12
    #sp 8->13  
    ind[1,0,0]=13    
    ind[2,0,0]=13   
    ind[3,0,0]=13                
    ind[4,0,0]=12            
    ind[5,0,0]=12            
    ind[6,0,0]=12                
    ind[7,0,0]=12            
    ind[8,0,0]=12                
    ind[4,1,0:2]=[10,11]            
    ind[5,1,0:2]=[10,11]            
    ind[6,1,0:2]=[10,11]            
    ind[7,1,0:2]=[10,11]            
    ind[8,1,0:2]=[10,11]            
    ind[4,2,0:2]=[10,11]            
    ind[5,2,0:2]=[10,11]            
    ind[6,2,0:2]=[10,11]            
    ind[7,2,0:2]=[10,11]            
    ind[8,2,0:2]=[10,11]            
    ind[4,3,0:2]=[10,11]            
    ind[5,3,0:2]=[10,11]            
    ind[6,3,0:2]=[10,11]            
    ind[7,3,0:2]=[10,11]            
    ind[8,3,0:2]=[10,11]        
        
    for i in range(noi):
        for j in range(noj):
            ht[i,j]=sum( [mat[i,j,k]*h[ind[i,j,k]] for k in range(cnt[i,j])] )
            st[i,j]=sum( [mat[i,j,k]*s[ind[i,j,k]] for k in range(cnt[i,j])] )     
            for a in range(3):
                dht[i,j,a]=sum( [mat[i,j,k]*dh[ind[i,j,k],a]+der[i,j,k,a]*h[ind[i,j,k]] for k in range(cnt[i,j])] )
                dst[i,j,a]=sum( [mat[i,j,k]*ds[ind[i,j,k],a]+der[i,j,k,a]*s[ind[i,j,k]] for k in range(cnt[i,j])] )          
    return ht, st, dht, dst                    


