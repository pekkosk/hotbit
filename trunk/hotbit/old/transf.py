from box import mix
s1='Au'
s2='C'
#dds ddp ddd pds pdp pps ppp sds sps sss
integrals0 =['dds','ddp','ddd','pds','pdp','pps','ppp','sds','sps','sss']
integrals1 =[                  'dps','dpp',            'dss','pss'      ]
angm=       [                  (2,1),(2,1),            (2,0),(1,0)      ]
ind1=       [                    3,    4,                7,    8        ]


t12=mix.find_value('%s_%s.par' %(s1,s2),'%s_%s_table' %(s1,s2),fmt='matrix')
t21=mix.find_value('%s_%s.par' %(s1,s2),'%s_%s_table' %(s2,s1),fmt='matrix')

r=t12[:,0]
N=len(r)

f=open('%s_%s_transformed.par' %(s1,s2),'w')
for H in [True,False]:
    add=[10,0][H]
    table, names=[], []
    
    # s1_s2 table as it is (but take only non-zero columns)
    for i0,name0 in enumerate(integrals0):
        ind=i0+1+add
        if max(abs(t12[:,ind]))<1E-10: continue
        table.append( t12[:,ind] )
        names.append( name0 )
            
    # add 4 additional integrals, if nonzero
    for i1 in range(4):
        ind=ind1[i1]+1+add
        if max(abs(t21[:,ind]))<1E-10: continue
        l1,l2=angm[i1]
        table.append( t21[:,ind]*(-1)**(l1+l2) )   
        names.append( integrals1[i1] )
            
    # write which tables are shown                           
    if H: print>>f, '%s_%s_H=\n            ' %(s1,s2),
    else: print>>f, '%s_%s_S=\n            ' %(s1,s2),
    for name in names:
        print>>f, '%14s' %name,
    print>>f
    for i in range(N):
        print>>f, '%14.6e' %r[i],
        for tab in table:
            print>>f, '%14.6e' %tab[i],
        print>>f
    print>>f, '\n\n\n'        
f.close()     
    
    
    

