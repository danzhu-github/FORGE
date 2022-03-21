import numpy as np
from scipy import stats

def pop(restart,nlons,nlats,ntime,eG_cst,eH_cst,c_dep,q,t2m,food_veg_ini,popu,fat,influx_fruit,decay_veg,food_veg,food_ani,fat_ani,tforage,eG,eH,tH,expend,\
estab_date,estab_elapse,estab_humsum,starv_humsum,fat_max_rec,fat_max_date):

  Nhum=1 # just one human funtional type
  estab_hum_a = 0.1
  estab_hum_b = 15.
  estab_hum_x0= 0.5
  expend_cst=8.37 #MJ/ind/day
  expend_veg=1.25 #MJ/ind/hr
  expend_ani=1.25 #MJ/ind/hr
  NE_ani = 9.8 #MJ/kgDM
  NE_veg = 9.8  #MJ/kgDM
  m_anabolism=54.6   # MJ(NE)/kg(fat) 
  m_catabolism=39.3  # MJ(NE)/kg(fat) 
  cdfnor_x_coe=0.
  cdfnor_sd_coe=0.125 
  hum_bw=50.  # body weight (kg/ind) 
  popu_init=1.e-9
  tforage_max=8. # hour/day
  tforage_min=0.5 # hour/day
  A=4000 # m2/hr/ind.

  #==============================================================
  #====== calculate some parameters =============================
  mort_hum = 1./80.   
  fat_max = hum_bw * 0.3
  cdfnor_x= cdfnor_x_coe * fat_max  
  cdfnor_sd= cdfnor_sd_coe * fat_max
  
  #==============================================================
  #====== initialization at beginning of run ====================
  valid=(t2m==t2m) # only land pixel
  
  estab_actual = np.zeros((nlats,nlons))*np.nan
  
  if restart=='N':
    popu[ valid] = popu_init
    fat[  valid] = 0.
    food_veg = food_veg_ini[0]*1
    tforage[valid]= tforage_max
    tH[valid] = tforage_max - tforage_min
    eG[valid] = eG_cst
    eH[valid] = eH_cst
    expend[valid]= expend_cst
    
    estab_date[valid]=364
    estab_elapse[valid]=0
    estab_humsum[valid] = 0
    starv_humsum[valid] = 0
    fat_max_rec[valid]  = -9999
    fat_max_date[valid] = 0
  
  popu_all    =np.zeros((ntime,Nhum,nlats,nlons))*np.nan 
  fat_all     =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  expend_all      =np.zeros((ntime,Nhum,nlats,nlons))*np.nan 
  estab_hum_all   =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  starv_hum_all   =np.zeros((ntime,Nhum,nlats,nlons))*np.nan 
  estab_actual_all=np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  intake_veg_all  =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  intake_ani_all  =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  mort_starv_all  =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  food_veg_all    =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  food_ani_all    =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  tforage_all     =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  tG_all          =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  tH_all          =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  eG_all          =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  eH_all          =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  benifit_ani_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  benifit_veg_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  
  estab_date_all   =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  estab_elapse_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  estab_humsum_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  starv_humsum_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  fat_max_rec_all  =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  fat_max_date_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan

  for itime in range(ntime):
  
      DoY= np.mod(itime,365)

      food_veg_all[    itime,0]=food_veg*1
      fat_old=fat*1

      tG=tforage-tH
      tforage0=tforage*1
      if np.any((tforage<=0)+(tG<=0)):
        print ('ERROR!!!!!!!!! tforage or tG<=0, itime=',itime,'MIN(tforage)=',np.nanmin(tforage),'loc=',np.nanargmin(tforage),'MIN(tG)=',np.nanmin(tG),'loc=',np.nanargmin(tG))
        sys.exit()
      if np.sum(valid)!=np.sum(tH==tH):
        print ('ERROR!!!!!!!!! (1) np.sum(valid)!=np.sum(tH==tH)',np.sum(valid),np.sum(tG==tG),np.sum(tH==tH),np.sum(tforage==tforage))
        sys.exit()
  
      # update foraging time
      if itime>0:
        cond_lazy = (fat>fat_max*0.5)*(delta_fat>0)
        cond_work=np.logical_not(cond_lazy)
        tforage[cond_lazy]=tforage[cond_lazy]-1
        tforage[cond_work]=tforage[cond_work]+1      

      tforage[tforage<tforage_min]=tforage_min
      tforage[tforage>tforage_max]=tforage_max

      tG=tG/tforage0*tforage
      tH=tH/tforage0*tforage

      # energy expenditure of foraging
      expend=expend_cst+ tG*expend_veg + tH*expend_ani

      benifit_veg=A * food_veg *eG * NE_veg
      benifit_ani=A * food_ani *eH * NE_ani

      # calculate daily intake: note that whenever fat is about to exceed fat_max, daily intake is reduced so that fat stays at fat_max
      tmp1= tG* benifit_veg
      tmp2= tH* benifit_ani

      tmp3 = tmp1+tmp2 # MJ/day/ind.

      cond1 =((fat+(tmp3-expend)/m_anabolism) >= fat_max) # approaching fat_max
      tmp3[cond1] =( (fat_max-fat[cond1])*m_anabolism+expend[cond1] ) # cannot exceed fat_max

      tmp4 = tmp3 * tmp1/(tmp1+tmp2) / NE_veg
      tmp5 = tmp3 * tmp2/(tmp1+tmp2) / NE_ani
      tmp4[(tmp1+tmp2)==0]=0.
      tmp5[(tmp1+tmp2)==0]=0.

      intake_veg =np.where(valid,tmp4,np.nan) # daily intake of plant food
      intake_ani =np.where(valid,tmp5,np.nan) # daily intake of animal food

      # reduce tG and tH in case humans are very full (fat is about to exceed fat_max)
      tG=np.where((food_veg>0.)*(intake_veg>1.e-6),intake_veg/food_veg/A/eG,tG)
      tH=np.where((food_ani>0.)*(intake_ani>1.e-6),intake_ani/food_ani/A/eH,tH)
      tG[tG<tforage_min]=tforage_min/2.
      tH[tH<tforage_min]=tforage_min/2.


      # intake cannot make food density below zero
      cond2=(food_veg - intake_veg*popu < 0.)
      tmp= food_veg / popu
      intake_veg[cond2]= tmp[cond2]*1

      cond3=(food_ani - intake_ani*popu < 0.)
      tmp= food_ani / popu
      intake_ani[cond3]= tmp[cond3]*1

      if np.any(intake_veg<0):
        print ('ERROR!!!!!!!!! intake_veg<0, itime=',itime,'MIN(intake_veg)=',np.nanmin(intake_veg),'loc=',np.nanargmin(intake_veg))
        sys.exit()
      if np.any(intake_ani<0):
        print ('ERROR!!!!!!!!! intake_ani<0, itime=',itime,'MIN(intake_ani)=',np.nanmin(intake_ani),'loc=',np.nanargmin(intake_ani))
        sys.exit()

      intake_total=intake_veg+intake_ani  # daily intake in dry mass unit
      intake_total_en=intake_veg*NE_veg+intake_ani*NE_ani  # daily intake in energy unit

      # update tforage, and tG/tH
      tmp=np.where(intake_total>0,intake_ani/intake_total,np.nan)
      meat_desire= np.exp((0.-q)*(tmp-1.)) 

      tforage= tG + tH
      if np.any((tforage<tforage_min/1.0001)+(tforage>tforage_max*1.0001)):
        print ('ERROR!!!!!!!!! tforage<tforage_min or >max, itime=',itime,'MIN(tforage)=',np.nanmin(tforage),'loc=',np.nanargmin(tforage),'MAX(tforage)=',np.nanmax(tforage),'loc=',np.nanargmax(tforage))
        sys.exit()

      tH_frac = tH/(tH+tG)

      tmp = benifit_ani + benifit_veg
      cond = (tmp>0)*(intake_total>0)
      tH_frac[cond] = benifit_ani[cond] / tmp[cond] * meat_desire[cond]

      tmp=np.where(intake_total>0,intake_ani/intake_total,np.nan)
      tH_frac=np.where(tmp<0.1,0.95,tH_frac)  # when meat fraction of diet is below 10%, tH_frac is set to 95%

      tH_frac[tH_frac<0.05]=0.05
      tH_frac[tH_frac>0.95]=0.95

      tH = tH_frac * tforage
      tG = (1.-tH_frac) * tforage

      if not np.allclose(tforage[valid],tG[valid]+tH[valid]):
        tmp=tforage-tG-tH; i=np.nanargmin(tmp); j=np.nanargmax(tmp)
        print ('ERROR!!!!!!!!! tG+tH!=tforage, itime=',itime,np.sum(valid),np.sum(tG==tG),np.sum(tH==tH),np.sum(tforage==tforage))
        sys.exit()

      # update eG and eH
      tmp = food_veg/(popu*50.)
      eG[valid]= eG_cst
      cond=(tmp==tmp)*(c_dep>0)
      eG[cond] = eG_cst*tmp[cond]/(tmp[cond]+c_dep)

      tmp = food_ani/(popu*50.)
      eH[valid]= eH_cst
      cond=(tmp==tmp)*(c_dep>0)
      eH[cond] = eH_cst*tmp[cond]/(tmp[cond]+c_dep)

      # update plant food density
      food_veg = food_veg - intake_veg*popu + influx_fruit[itime]
      food_veg = food_veg - decay_veg[itime]*food_veg
      if np.any(food_veg<0):
        print ('ERROR!!!!!!!!! food_veg<0, itime=',itime,'MIN(food_veg)=',np.nanmin(food_veg),'loc=',np.nanargmin(food_veg))
        sys.exit()

      # update animal food density
      food_ani = food_ani - intake_ani*popu
      if np.any(food_ani<0):
        print ('ERROR!!!!!!!!! food_ani<0, itime=',itime,'MIN=',np.nanmin(food_ani),'loc=',np.nanargmin(food_ani))
        sys.exit()

      # update fat reserve
      m_use=np.where(intake_total_en >= expend,m_anabolism,m_catabolism)
      fat = fat + (intake_total_en - expend)/m_use
      if np.any(fat>1.0001*fat_max):
        print ('ERROR!!!!!!!! fat>fat_max, itime=',itime,'MAX(fat)=',np.nanmax(fat),'loc=',np.nanargmax(fat))
        sys.exit()

      # record the date with maximum body fat, on which the birth will take place
      tmp=(fat>=fat_max_rec)
      fat_max_date[tmp]=DoY
      fat_max_rec[tmp] = fat[tmp]*1

      # calculate birth and mortality rate  
      estab_hum =estab_hum_a/(1.+np.exp(-estab_hum_b*(fat/fat_max-estab_hum_x0)))
  
      starv_hum = np.zeros(fat.shape)*np.nan
      for i in range(nlats):
          for j in range(nlons):
                  starv_hum[i,j]=stats.norm.cdf(cdfnor_x,fat[i,j],cdfnor_sd)
  
      # temporary variables to calcuate annual mean birth and mortality rate later
      estab_humsum=estab_humsum+estab_hum   
      starv_humsum=starv_humsum+starv_hum

      # START birth date
      e=(DoY==estab_date)
      if np.any(estab_elapse[e]==0):
        print ('ERROR!!!!!!!!!!! check estab_elapse')
        sys.exit()
      e=(DoY==estab_date)*(estab_elapse>330)

      estab_actual[e]=estab_humsum[e]/estab_elapse[e]
      fat[e]=fat[e]-estab_actual[e]*(313.*4.184/1000.*estab_elapse[e])/m_catabolism

      tmp=(estab_actual<0)
      if np.any(tmp[e]==1):
        print ('ERROR!!!!!!!!!!! estab_actual is negative')
        sys.exit()

      estab_actual_all[itime,0][e]=estab_actual[e]*1
      mort_starv_all[  itime,0][e]=starv_humsum[e]/estab_elapse[e]
 
      cond_low=(estab_actual-mort_hum-starv_humsum <= 0.)*e 
      ncond=   (estab_actual-mort_hum-starv_humsum >  0.)*e
      popu[ncond]=popu[ncond]*(1+estab_actual[ncond]-mort_hum-starv_humsum[ncond]/estab_elapse[ncond])
      popu[cond_low]=popu[cond_low]*(1+estab_actual[cond_low]-starv_humsum[cond_low]/estab_elapse[cond_low])
  
      estab_elapse[e]=0
      estab_humsum[e]=0
      starv_humsum[e]=0

      estab_date[e]=fat_max_date[e]*1
      fat_max_rec[e]=fat[e]*1

      # END birth date

      # reset to initial if popu too low
      reset=(popu<popu_init)
      fat[reset*(fat<0)]=0.
      popu[reset]=popu_init

      estab_elapse=estab_elapse+1
 
      delta_fat=fat - fat_old 
  
      popu_all[    itime,0]=popu    *1
      fat_all[     itime,0]=fat     *1  
      expend_all[      itime,0]=expend      *1 
      estab_hum_all[   itime,0]=estab_hum   *1    
      starv_hum_all[   itime,0]=starv_hum   *1 
      intake_veg_all[  itime,0]=intake_veg*1
      intake_ani_all[  itime,0]=intake_ani*1
      tforage_all[     itime,0]=tforage *1  
      tG_all[          itime,0]=tG *1  
      tH_all[          itime,0]=tH *1  
      eG_all[          itime,0]=eG *1  
      eH_all[          itime,0]=eH *1  
      food_ani_all[    itime,0]=food_ani *1
      benifit_veg_all[ itime,0]=benifit_veg *1  
      benifit_ani_all[ itime,0]=benifit_ani *1  

      estab_date_all[  itime,0]=estab_date  *1
      estab_elapse_all[itime,0]=estab_elapse*1
      estab_humsum_all[itime,0]=estab_humsum*1
      starv_humsum_all[itime,0]=starv_humsum*1
      fat_max_date_all[itime,0]=fat_max_date*1
      fat_max_rec_all[ itime,0]=fat_max_rec *1

      if itime == ntime-1:  food_veg_all[    itime,0]=food_veg*1  # this is to avoid drift in food_veg,if one-year forcing is looped

  ###################################################################
  return popu_all,fat_all,expend_all,intake_veg_all,estab_hum_all,starv_hum_all,estab_actual_all,mort_starv_all,food_veg_all,tforage_all,tG_all,eG_all,eH_all,intake_ani_all,food_ani_all,tH_all,benifit_veg_all,benifit_ani_all,estab_date_all,estab_elapse_all,estab_humsum_all,starv_humsum_all,fat_max_rec_all,fat_max_date_all
