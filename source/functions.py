import numpy as np


def time_aggre(res,varin):

    nds=varin.shape
    nyr=int(nds[0]/365)

    if res=='daily': varout = varin * 1
    if res=='month':
        a=varin * 1
        a.shape=(nyr,365,nds[1],nds[2],nds[3])
        b=np.zeros((nyr,12,nds[1],nds[2],nds[3]))*np.nan    
        ind=np.cumsum(calendar.mdays)  
        for imon in range(12):   
            b[:,imon]=np.nanmean(a[:,ind[imon]:ind[imon+1]],axis=1)
        varout=b.reshape(nyr*12,nds[1],nds[2],nds[3])
    if res=='year':
        a=varin * 1
        a.shape=(nyr,365,nds[1],nds[2],nds[3])
        varout=np.nanmean(a,axis=1)

    return varout
