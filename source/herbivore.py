import numpy as np
from scipy import stats

def animal(restart,nlons,nlats,ntime,Nani,t2m,bm_avail_ini,lit_avail_ini,influx_bm,decay_bm,influx_lit,decay_lit,bm_avail,lit_avail,\
estab_date,estab_elapse,estab_anisum,starv_anisum,fat_max_rec,fat_max_date,\
ani_popu,fat_ani):

  dt = 1 
  estab_ani_a =0.2
  estab_ani_b =10.
  estab_ani_x0=0.5 
  expend_a=2. 
  expend_b=0.0079 
  ne_AGB=5.46      # MJ/kgDM 
  ne_litter=3.12   # MJ/kgDM 
  m_anabolism=54.6   # MJ(NE)/kg(fat) 
  m_catabolism=39.3  # MJ(NE)/kg(fat) 
  delai_ugb_max=-5. 
  cdfnor_x_coe=-0.2
  cdfnor_sd_coe=0.125 
  t2m_lt_min=-5. 
  ani_bw=180.  # body weight (kg/ind) 
  KK=18000./ani_bw
  ani_popu_init=1.e-9
  Ldecay = 0.96 # decay rate of the edible litter

  #==============================================================
  #====== calculate some parameters =============================
  intake_AGB_max =    0.034*np.exp(3.57*0.7)*ani_bw**(0.077*np.exp(0.7)+0.73)   # MJ/day/head
  intake_litter_max = 0.034*np.exp(3.57*0.4)*ani_bw**(0.077*np.exp(0.4)+0.73)   # MJ/day/head
  try_able_grazing = intake_AGB_max/ne_AGB+1 
  mort_ani = 1./25.   
  fat_ani_max = ani_bw * 0.3 
  cdfnor_x= cdfnor_x_coe * fat_ani_max  
  cdfnor_sd= cdfnor_sd_coe * fat_ani_max
  
  #==============================================================
  #====== initialization at beginning of run ====================
  valid=(t2m==t2m) # only land pixel
  
  delai_ugb    = np.zeros((Nani,nlats,nlons))*np.nan
  able_grazing = np.zeros((Nani,nlats,nlons))*np.nan
  ugb          = np.zeros((Nani,nlats,nlons))*np.nan
  estab_actual = np.zeros((Nani,nlats,nlons))*np.nan
  delai_ugb[valid] = delai_ugb_max

  if restart=='N':
    ani_popu[ valid] = ani_popu_init
    fat_ani[  valid] = 0.
    bm_avail = bm_avail_ini*1
    lit_avail =lit_avail_ini*1
    estab_date[valid]=364
    estab_elapse[valid]=0
    estab_anisum[valid] = 0
    starv_anisum[valid] = 0
    fat_max_rec[valid]  = -9999
    fat_max_date[valid] = 0
  
  ani_popu_all    =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  fat_ani_all     =np.zeros((ntime,Nani,nlats,nlons))*np.nan   
  delai_ugb_all   =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  ugb_all         =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  expend_all      =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  estab_ani_all   =np.zeros((ntime,Nani,nlats,nlons))*np.nan   
  starv_ani_all   =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  estab_actual_all=np.zeros((ntime,Nani,nlats,nlons))*np.nan
  intake_ind_all  =np.zeros((ntime,Nani,nlats,nlons))*np.nan   
  intake_ind_litter_all  =np.zeros((ntime,Nani,nlats,nlons))*np.nan   
  bm_avail_all    =np.zeros((ntime,Nani,nlats,nlons))*np.nan
  lit_avail_all   =np.zeros((ntime,Nani,nlats,nlons))*np.nan
  estab_date_all  =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  estab_elapse_all =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  estab_anisum_all =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  starv_anisum_all =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  fat_max_rec_all  =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  fat_max_date_all =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  
  for itime in range(ntime):

      DoY= np.mod(itime,365)

      bm_avail_all[    itime]=bm_avail*1
      lit_avail_all[   itime]=lit_avail*1
  
      expend=np.where(t2m<(t2m_lt_min+273.15),\
                 expend_a/np.exp(expend_b*t2m_lt_min)*(ani_bw*1000.)**0.75/1000.,\
                 expend_a/np.exp(expend_b*(t2m-273.15))*(ani_bw*1000.)**0.75/1000.)
      expend[np.logical_not(valid)]=np.nan

      able_grazing = ani_popu*try_able_grazing

      cond1=(bm_avail>=able_grazing) # enough fresh biomass
      delai_ugb[valid*cond1]=delai_ugb[valid*cond1]+1
      delai_ugb[valid*(np.logical_not(cond1))]=delai_ugb_max
  
      cond2=(delai_ugb>=0) # after a sufficient time (5 days)
      cond3=(lit_avail>=able_grazing) # enough dead biomass
  
      # determine what to eat
      ugb[valid]=0
      ugb[cond1*cond2]=2
      ugb[(ugb!=2)*cond3]=1
  
      # start the eating
      intake_ind       =np.where(valid,0.,np.nan)
      intake_ind_litter=np.where(valid,0.,np.nan)
  
      cond1=(ugb==2)*(ani_popu>0)
      cond2=(ugb==1)*(ani_popu>0)
      ncond=np.logical_not(cond1+cond2)*valid
  
      intake_ind[cond1] = intake_AGB_max/ne_AGB
      cond3 =((fat_ani+(intake_ind*ne_AGB-expend)/m_anabolism) >= fat_ani_max) # approaching fat_ani_max
      intake_ind[cond3] =( (fat_ani_max-fat_ani[cond3])*m_anabolism+expend[cond3] )/ne_AGB # cannot exceed fat_ani_max
      if np.any(intake_ind<0):
        print ('ERROR! intake_ind<0, itime=',itime)
        sys.exit()
  
      intake_ind_litter[cond2] = intake_litter_max/ne_litter

      # update biomass and litter pool
      bm_avail = bm_avail - intake_ind*ani_popu + influx_bm[itime]
      bm_avail = bm_avail - decay_bm[itime]*bm_avail
      if np.any(bm_avail<0):
        print ('ERROR!!!!!!!!! bm_avail<0, itime=',itime,'MIN(bm_avail)=',np.nanmin(bm_avail),'loc=',np.nanargmin(bm_avail))
        sys.exit()
  
      lit_avail = lit_avail - intake_ind_litter*ani_popu + influx_lit[itime]
      lit_avail = lit_avail - decay_lit[itime]*lit_avail
      lit_avail = Ldecay*lit_avail
      if np.any(lit_avail<0):
        print ('ERROR!!!!!!!!! lit_avail<0, itime=',itime,'MIN(lit_avail)=',np.nanmin(lit_avail),'loc=',np.nanargmin(lit_avail))
        sys.exit()
  
      # update state variables
      fat_ani[cond1]=fat_ani[cond1] + (intake_ind[cond1]*ne_AGB-expend[cond1])/m_anabolism
      fat_ani[cond2]=fat_ani[cond2] + (intake_ind_litter[cond2]*ne_litter-expend[cond2])/m_catabolism
      fat_ani[ncond]=fat_ani[ncond]-expend[ncond]/m_catabolism
      if np.any(fat_ani>fat_ani_max):
        print ('ERROR! fat_ani>fat_ani_max, itime=',itime)
        sys.exit()

      # record the maximum fat date
      tmp=(fat_ani>=fat_max_rec)
      fat_max_date[tmp]=DoY
      fat_max_rec[tmp] = fat_ani[tmp]*1
 
      # calculate birth rate and mortality rate 
      estab_ani =estab_ani_a/(1.+np.exp(-estab_ani_b*(fat_ani/fat_ani_max-estab_ani_x0)))
      estab_ani[fat_ani<=0]=0.
  
      starv_ani =np.where(fat_ani>=0,0.,np.nan)
      for k in range(Nani):
        for i in range(nlats):
            for j in range(nlons):
                if fat_ani[k,i,j]<0:
                    starv_ani[k,i,j]=stats.norm.cdf(cdfnor_x,fat_ani[k,i,j],cdfnor_sd)
  
      estab_anisum=estab_anisum+estab_ani      
      starv_anisum=starv_anisum+starv_ani
  
      # START birth date
      e=(DoY==estab_date)
      if np.any(estab_elapse[e]==0):
        print ('ERROR!!!!!!!!!!! check estab_elapse')
        sys.exit()
      e=(DoY==estab_date)*(estab_elapse>330)

      cond1=(fat_ani*m_catabolism/expend > estab_anisum )*e
      cond2=(fat_ani*m_catabolism/expend <=estab_anisum )*e

      estab_actual[cond1]=estab_anisum[cond1]/estab_elapse[cond1]
      estab_actual[cond2]=fat_ani[cond2]*m_catabolism/(expend[cond2]*estab_elapse[cond2])
      estab_actual[estab_actual<0]=0.
  
      fat_ani[e]=fat_ani[e]-estab_actual[e]*(expend[e]*estab_elapse[e])/m_catabolism
  
      ani_popu[e]=ani_popu[e]*(1+estab_actual[e]-mort_ani-starv_anisum[e]/estab_elapse[e])-(estab_ani_a-mort_ani)/KK*ani_popu[e]*ani_popu[e]*1.e6

      estab_elapse[e]=0
      estab_anisum[e]=0
      starv_anisum[e]=0

      estab_date[e]=fat_max_date[e]*1
      fat_max_rec[e]=fat_ani[e]*1

      # END birth date

      # reset to initial if popu too low
      reset=(ani_popu<ani_popu_init)  
      fat_ani[reset]=0.
      ani_popu[reset]=ani_popu_init

      estab_elapse=estab_elapse+1  
  
      ani_popu_all[    itime]=ani_popu    *1
      fat_ani_all[     itime]=fat_ani     *1  
      delai_ugb_all[   itime]=delai_ugb   *1 
      ugb_all[         itime]=ugb         *1 
      expend_all[      itime]=expend      *1 
      estab_ani_all[   itime]=estab_ani   *1    
      starv_ani_all[   itime]=starv_ani   *1 
      estab_actual_all[itime]=estab_actual*1
      intake_ind_all[  itime]=intake_ind*1
      intake_ind_litter_all[itime]=intake_ind_litter*1
      estab_date_all[  itime]=estab_date  *1
      estab_elapse_all[itime]=estab_elapse*1
      estab_anisum_all[itime]=estab_anisum*1
      starv_anisum_all[itime]=starv_anisum*1
      fat_max_date_all[itime]=fat_max_date*1
      fat_max_rec_all[ itime]=fat_max_rec *1

  ###################################################################
  return ani_popu_all,fat_ani_all,delai_ugb_all,ugb_all,expend_all,intake_ind_all,intake_ind_litter_all,estab_ani_all,starv_ani_all,estab_actual_all,bm_avail_all,lit_avail_all,estab_date_all,estab_elapse_all,estab_anisum_all,starv_anisum_all,fat_max_rec_all,fat_max_date_all
