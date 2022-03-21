import numpy as np
from netCDF4 import Dataset

def write_var(res,pathin,pathout,filein,fileout,nlat,nlon,ntime,Nveg,Nani,Nhum,outlist,varout):
  
  ncref=Dataset(pathin+filein,'r')
  nc=Dataset(pathout+fileout,'w',format='NETCDF3_CLASSIC')
  
  nc.createDimension('lat',nlat)
  nc.createDimension('lon',nlon)
  nc.createDimension('veg',Nveg)
  nc.createDimension('ani',Nani)
  nc.createDimension('hum',Nhum)
  nc.createDimension('time',None)
  
  var=nc.createVariable('lat','f',('lat',),zlib=True)
  var[:]=ncref.variables['lat'][:]
  vattr=['standard_name','long_name','units','axis']
  for a in vattr:
    attr=ncref.variables['lat'].__getattribute__(a)
    nc.variables['lat'].__setattr__(a,attr)
  
  var=nc.createVariable('lon','f',('lon',),zlib=True)
  var[:]=ncref.variables['lon'][:]
  vattr=['standard_name','long_name','units','axis']
  for a in vattr:
    attr=ncref.variables['lon'].__getattribute__(a)
    nc.variables['lon'].__setattr__(a,attr)
  
  var=nc.createVariable('time','f',('time',),zlib=True)
  if res=='daily': var[:]=np.arange(ntime)*86400
  if res=='month':
      tmp=np.cumsum(calendar.mdays)
      var[:]=np.array([ (i/12)*86400*365+(tmp[np.mod(i,12)]+15)*86400 for i in range(ntime)])
  if res=='year' : var[:]=np.arange(ntime)*86400*365+86400*182
  nc.variables['time'].__setattr__('units',u'seconds since 1901-01-01 00:00:00')
  nc.variables['time'].__setattr__('calendar',u'365_day')
  nc.variables['time'].__setattr__('long_name',u'Time axis')
  nc.variables['time'].__setattr__('standard_name',u'time')
  nc.variables['time'].__setattr__('axis',u'T')
  
  #####################################################################
  for loop in range(len(outlist)):
      tmp=varout[outlist[loop][0]]*1
      tmp2=outlist[loop][1]
      name=outlist[loop][0]
      if tmp2=='2':    var=nc.createVariable(name,'f',('lat','lon'),zlib=True)
      if tmp2=='3':    var=nc.createVariable(name,'f',('time','lat','lon'),zlib=True)
      if tmp2=='veg':  var=nc.createVariable(name,'f',('time','veg','lat','lon'),zlib=True)
      if tmp2=='ani':  var=nc.createVariable(name,'f',('time','ani','lat','lon'),zlib=True)
      if tmp2=='hum':  var=nc.createVariable(name,'f',('time','hum','lat','lon'),zlib=True)

      tmp[tmp!=tmp]=-9.e20
      var[:]=tmp

      nc.variables[name].__setattr__('missing_value',-9.e20)
  #####################################################################
  ncref.close()
  nc.close()

  return
