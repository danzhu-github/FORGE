import numpy as np
import readNC
import writeNC
import herbivore 
import human
import functions as Fun


def run(pathin,pathout,res,run_year,phi_v,phi_a,eG_cst,eH_cst,c_dep,q,filein_pre,fileout):
  
  ### variable name in the input file
  varlist=['lon','lat','veg','time',\
          't2m','prec', 'land_cov','influx_bm_Nani','influx_lit_Nani','decay_bm_Nani','decay_lit_Nani','npp',      'decay_bm_Nhum']
  ### variable name used in the FORGE model,corresponding to the names from the input
  varuse =['',    '',   '' ,  '',           \
          't2m','prec', 'land_cov','in_bm_ani',     'in_lit_ani',     'turn_bm_ani',  'turn_lit_ani',  'in_bm_hum','turn_bm_hum']

  Nani=2 # two herbivore funtional types (HFTs): grazers and browsers
  Nhum=1 # just one human funtional type
  CtoDM=0.45  # conversion factor between dry mass and carbon mass
  edi_frac_tree=0.1 # 10% of tree leaf and fruits (edible to herbivores) is reachable to browsers
  ani_edi=0.28 # conversion between live weight and dry mass after excluding water and bones for the animal food

  ### variable list to be outputted
  adict=locals()
  outlist=[['ani_popu','ani'],['fat_ani','ani'],['delai_ugb','ani'],['ugb_ani','ani'],['expend_ani','ani'],['intake_ind_ani','ani'],['intake_ind_litter_ani','ani'],['estab_ani','ani'],['starv_ani','ani'],['estab_actual_ani','ani'],['bm_avail','ani'],['lit_avail','ani'],['est_date','ani'],['est_elapse','ani'],['in_bm_avail','ani'],['in_lit_avail','ani'], \
           ['hum_popu','hum'],['fat_hum','hum'],['expend_hum','hum'],['intake_veg','hum'],['estab_hum','hum'],['starv_hum','hum'],['estab_actual_hum','hum'],['mort_starv_hum','hum'],['food_veg','hum'],\
           ['tforage','hum'],['tG','hum'],['eG','hum'],['eH','hum'],['intake_ani','hum'],['food_ani','hum'],['tH','hum'],['benifit_veg','hum'],['benifit_ani','hum'],['est_date_hum','hum'],['est_elapse_hum','hum'],\
           ['influx_v','hum']]
  
  DayInYr = 365 # days in a year
  ntime_run = 365*run_year
  
  for loop in range(run_year):
  
    print (loop)
    #==============================================================
    #=== read inputs from netcdf files ============================
    #==============================================================
    input_year = 2001+np.mod(loop,10)  # the 10-yr inputs files are cycled
    filein = filein_pre+str(input_year)+'.nc'
  
    nlons,nlats,nlevs,ntime,varin = readNC.read_var(pathin,filein,varlist,varuse)
    if loop==0: print ('nlats=',nlats,'nlons=',nlons,'ntime=',ntime,'ntime_run=',ntime_run)
    
    if np.mod(ntime_run,ntime)!=0 or DayInYr!=ntime: 
        print ('WRONG in ntime_run!')
        sys.exit()
    
    #==============================================================
    #=== prepare inputs for the herbiovore module =================
    #==============================================================
    fgrass= varin['land_cov'][:,0]*1  # grass fractional cover
    ftree = varin['land_cov'][:,1]*1  # tree fractional cover
    fgrass_yr= np.mean(fgrass,axis=0) # annual mean
    ftree_yr = np.mean(ftree, axis=0) # annual mean
  
    t2m_use=np.repeat(varin['t2m'].reshape(1,nlats*nlons),[Nani]*1,axis=0)  # to add the HFT dimension 
    t2m_use.shape=(Nani,nlats,nlons)
  
    bm_avail_ini_Nani = np.zeros((Nani,nlats,nlons))  # initial food (plant living biomass) density for herbivores (set to zero)
    lit_avail_ini_Nani = np.zeros((Nani,nlats,nlons)) # initial food (plant litter) density for herbivores (set to zero)
  
    influx_bm_avail_Nani= varin['in_bm_ani']/(1000. * CtoDM) # influx to the edible biomass for herbivores, unit: kgDM/m2/day
    influx_bm_avail_Nani[:,1,:,:] = influx_bm_avail_Nani[:,1,:,:] *edi_frac_tree  # for browsers, only 10% is reachable.

    decay_bm_avail_Nani= varin['turn_bm_ani']*1  # daily decay rate of the edible biomass for herbivores
  
    influx_lit_avail_Nani= varin['in_lit_ani']/(1000. * CtoDM) # influx to the edible litter for herbivores, unit: kgDM/m2/day
    influx_lit_avail_Nani[:,1,:,:] = influx_lit_avail_Nani[:,1,:,:] *edi_frac_tree
  
    decay_lit_avail_Nani= varin['turn_lit_ani']*1 # daily decay rate of the edible litter for herbivores

    ### do some check ###
    if np.any(influx_bm_avail_Nani<0):
      loc=np.nanargmin(influx_bm_avail_Nani)
      print ('influx_bm_avail_Nani<0 MIN=',np.nanmin(influx_bm_avail_Nani),'loc=',loc)
      sys.exit()
    if np.any(influx_lit_avail_Nani<0):
      loc=np.nanargmin(influx_lit_avail_Nani)
      print ('influx_lit_avail_Nani<0 MIN=',np.nanmin(influx_lit_avail_Nani),'loc=',loc)
      sys.exit()
    if np.any((decay_bm_avail_Nani<0)+(decay_bm_avail_Nani>1)):
      print ('decay_bm_avail_Nani<0 or >1')
      sys.exit()
    if np.any((decay_lit_avail_Nani<0)+(decay_lit_avail_Nani>1)):
      print ('decay_lit_avail_Nani<0 or >1')
      sys.exit()

    #==============================================================
    #=== calculate inputs for the human module ====================
    #==============================================================
    food_veg_ini=np.zeros((nlats,nlons))
 
    ### calculate the influx of plant food for humans, unit: kgDM/m2/day. see Equation(7)
    npp = varin['in_bm_hum'][:,:,:,:]/(1000. * CtoDM)
    tmp=varin['prec'][:,:]*1
    fBNPP=(88.3-0.0534*tmp)/100. ; fBNPP[fBNPP<0.2]=0.2
    influx_veg = ( (npp[:,0]*fgrass + npp[:,1]*ftree)*0.065 + npp[:,0]*fgrass*fBNPP ) * phi_v

    decay_veg= varin['turn_bm_hum'][:,0,:,:]*1  # daily decay rate of the edible plants for humans
  
    ### do some check ###
    if np.any(influx_veg<0):
      loc=np.nanargmin(influx_veg)
      print ('influx_veg<0 MIN=',np.nanmin(influx_veg),'loc=',loc)
      sys.exit()
    if np.any((decay_veg<0)+(decay_veg>1)):
      print ('decay_veg<0 or >1')
      sys.exit()

    #==============================================================
    #=== run the herbivore and human model ========================
    #==============================================================
    if loop==0: # at the beginning of simulation, do the initialization 
      ani_popu_0= np.zeros((Nani,nlats,nlons))*np.nan
      fat_ani_0 = np.zeros((Nani,nlats,nlons))*np.nan
      bm_avail_0= np.zeros((Nani,nlats,nlons))*np.nan
      lit_avail_0=np.zeros((Nani,nlats,nlons))*np.nan
      estab_date_0=np.zeros((Nani,nlats,nlons))*np.nan
      estab_elapse_0=np.zeros((Nani,nlats,nlons))*np.nan
      estab_anisum_0=np.zeros((Nani,nlats,nlons))*np.nan
      starv_anisum_0=np.zeros((Nani,nlats,nlons))*np.nan
      fat_max_rec_0=np.zeros((Nani,nlats,nlons))*np.nan
      fat_max_date_0=np.zeros((Nani,nlats,nlons))*np.nan
      
      hum_popu_0= np.zeros((nlats,nlons))*np.nan
      fat_hum_0 = np.zeros((nlats,nlons))*np.nan
      food_veg_0= np.zeros((nlats,nlons))*np.nan
      tforage_0 = np.zeros((nlats,nlons))*np.nan
      eG_0      = np.zeros((nlats,nlons))*np.nan
      eH_0      = np.zeros((nlats,nlons))*np.nan
      tH_0      = np.zeros((nlats,nlons))*np.nan
      expend_0  = np.zeros((nlats,nlons))*np.nan
      estab_date_hum_0  =np.zeros((nlats,nlons))*np.nan
      estab_elapse_hum_0=np.zeros((nlats,nlons))*np.nan
      estab_humsum_0    =np.zeros((nlats,nlons))*np.nan
      starv_humsum_0    =np.zeros((nlats,nlons))*np.nan
      fat_max_rec_hum_0 =np.zeros((nlats,nlons))*np.nan
      fat_max_date_hum_0=np.zeros((nlats,nlons))*np.nan
  
      restart = 'N'
  
    else:
      restart = 'Y'  

    ### run the herbivore module
    Vs = herbivore.animal(restart,nlons,nlats,ntime,Nani,t2m_use,bm_avail_ini_Nani,lit_avail_ini_Nani,influx_bm_avail_Nani,decay_bm_avail_Nani,influx_lit_avail_Nani,decay_lit_avail_Nani,bm_avail_0,lit_avail_0,\
         estab_date_0,estab_elapse_0,estab_anisum_0,starv_anisum_0,fat_max_rec_0,fat_max_date_0,\
         ani_popu_0,fat_ani_0)
    ### record the state variables
    ani_popu_0=Vs[0][-1]*1
    fat_ani_0 =Vs[1][-1]*1
    bm_avail_0=Vs[10][-1]*1
    lit_avail_0=Vs[11][-1]*1
    estab_date_0   = Vs[12][-1]*1
    estab_elapse_0 = Vs[13][-1]*1
    estab_anisum_0 = Vs[14][-1]*1
    starv_anisum_0 = Vs[15][-1]*1
    fat_max_rec_0  = Vs[16][-1]*1
    fat_max_date_0 = Vs[17][-1]*1
   
    ### prepare the output variables
    s=influx_bm_avail_Nani.shape ; s=(1,s[0],s[1],s[2],s[3])
    Vs_out=np.concatenate((Vs[:14],influx_bm_avail_Nani.reshape(s),influx_lit_avail_Nani.reshape(s)),axis=0) # write out the food influx for herbivores as well

    ### convert from kg(live weight)/m2 to edible/accessible kgDM/m2 
    ani_tmp1= ani_popu_0[0] * fgrass_yr*180.*ani_edi * phi_a
    ani_tmp2= ani_popu_0[1] * ftree_yr *180.*ani_edi * phi_a
    ani_DM_0 = ani_tmp1+ani_tmp2 # sum up the density on tree and grass fractional covers to derive the density for the whole grid cell
    ani_old = ani_DM_0 *1
   
    ### run the human module
    Ps  = human.pop(restart,nlons,nlats,ntime,eG_cst,eH_cst,c_dep,q,varin['t2m'],food_veg_ini,hum_popu_0,fat_hum_0,influx_veg,decay_veg,food_veg_0,ani_DM_0,fat_ani_0,tforage_0,eG_0,eH_0,tH_0,expend_0,\
  estab_date_hum_0,estab_elapse_hum_0,estab_humsum_0,starv_humsum_0,fat_max_rec_hum_0,fat_max_date_hum_0)
    
    ### record the state variables
    hum_popu_0= Ps[0][-1,0]*1
    fat_hum_0 = Ps[1][-1,0]*1
    food_veg_0= Ps[8][-1,0]*1
    tforage_0 = Ps[9][-1,0]*1
    eG_0      =Ps[11][-1,0]*1
    eH_0      =Ps[12][-1,0]*1
    tH_0      =Ps[15][-1,0]*1
    ani_DM_0  =Ps[14][-1,0]*1
    estab_date_hum_0  = Ps[18][-1,0]*1
    estab_elapse_hum_0= Ps[19][-1,0]*1
    estab_humsum_0    = Ps[20][-1,0]*1
    starv_humsum_0    = Ps[21][-1,0]*1
    fat_max_rec_hum_0 = Ps[22][-1,0]*1
    fat_max_date_hum_0= Ps[23][-1,0]*1
  
    ### prepare the output variables
    s=influx_veg.shape ; s=(1,s[0],1,s[1],s[2])
    Ps_out=np.concatenate((Ps[:20],influx_veg.reshape(s)),axis=0) # write out the influx to plant food for humans as well
    
    ### to subtract the hunted herbivores: assume a proportional reduction between grazers and browsers
    ani_new = ani_DM_0 *1
    if np.any(ani_new>ani_old):
      print ('ERROR!!!!!!!!! animal increases after hunting ')
      sys.exit()
    ani_dif1=np.where(ani_old==0,0,ani_tmp1*(ani_old-ani_new)/ani_old)
    ani_dif2=np.where(ani_old==0,0,ani_tmp2*(ani_old-ani_new)/ani_old)
 
    ani_popu_0[0]=ani_popu_0[0] - np.where(fgrass_yr>0, ani_dif1/(fgrass_yr*180.*ani_edi), 0.) # note: the unit of animal population density in the herbivore module is per unit PFT (grass or tree) area, instead of per ground area, so here it needs to divide by the fractional cover of grass or tree before being used in the herbivore module
    ani_popu_0[1]=ani_popu_0[1] - np.where(ftree_yr >0, ani_dif2/(ftree_yr *180.*ani_edi), 0.)
    if np.any(ani_popu_0<0):
      print ('ERROR!!!!!!!!! animal less than 0')
      sys.exit()
    
    ### prepare the output list: average along the time dimension if res='month' or 'year'
    Numv=16
    if loop==0:
      for i in range(len(outlist)):
        if i<Numv:  
            adict[outlist[i][0] ]=Fun.time_aggre(res,Vs_out[i] )
        if i>=Numv: 
            adict[outlist[i][0] ]=Fun.time_aggre(res,Ps_out[i-Numv] )
    else:
      for i in range(len(outlist)):
        if i<Numv:  tmp=Fun.time_aggre(res,Vs_out[i] )
        if i>=Numv: tmp=Fun.time_aggre(res,Ps_out[i-Numv] )
        adict[outlist[i][0]] =np.concatenate((adict[outlist[i][0]],tmp),axis=0)
    
  
  #==============================================================
  #=== write outputs into netcdf files ==========================
  #==============================================================
  varout={outlist[0][0]: adict[outlist[0][0]]} # a dictionary to store all the variables  
  for i in range(1,len(outlist)): 
      varout.update({outlist[i][0]: adict[outlist[i][0]]}) #update the dictionary
  
  if res=='daily': ntime_out=ntime_run
  if res=='month': ntime_out=ntime_run/365*12
  if res=='year':  ntime_out=ntime_run/365
  
  writeNC.write_var(res,pathin,pathout,filein,fileout,nlats,nlons,ntime_out,nlevs,Nani,Nhum,outlist,varout)
  
