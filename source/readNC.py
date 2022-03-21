import numpy as np
from netCDF4 import Dataset

def read_var(path,filename,varlist,varuse):

    adict=locals()
    nc=Dataset(path+filename,'r')
    nlons=nc.dimensions[varlist[0]].size
    nlats=nc.dimensions[varlist[1]].size
    nlevs=nc.dimensions[varlist[2]].size
    ntime=nc.dimensions[varlist[3]].size

    for i in range(4,len(varlist)):
        adict[varuse[i]]=nc.variables[varlist[i]][:]
        adict[varuse[i]]=adict[varuse[i]].filled(np.nan)
        adict[varuse[i]][adict[varuse[i]]<-8e+20]=np.nan# NOTE
    nc.close()

    varin={varuse[4]: adict[varuse[4]]} # a dictionary to store all the variables
    for i in range(5,len(varlist)):
        varin.update({varuse[i]: adict[varuse[i]]}) #update the dictionary
    
    return nlons,nlats,nlevs,ntime,varin
