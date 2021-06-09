import matplotlib.pyplot as plt 
from astropy.table import Table, Column
from astropy.io import fits
import getweather as gw
import getrank as gr
import telescope as ts
from sitesky import *
import numpy as np
import logging

# revised in Aug 31 2020, adding additional outputs
#          such as airmass, v-band sky brightness and seeing
# read the fields
hdu = fits.open('datacube.fits')
data = hdu[1].data
ras = np.array(data['raa'])
decs = np.array(data['deca'])
fieldid = np.arange(0,len(ras))
nfield = len(ras)

# read the sun set and sun rise time
hdu = fits.open('sunsetrise'+'%4i' % conf['year'] +'.fits')
tsunset = hdu[1].data
tseta = np.array(tsunset['tset'])
trisa = np.array(tsunset['trise'])

# read the avaiable time
tmp = fits.open('fordataresult.fits')
t1s = tmp[1].data
avtaa = t1s['avati']

# read the weather data
filename='lijiangmodel.fits'
tmp = fits.open(filename)
tweather = tmp[1].data

if isLeap(conf['year']):
    nial = 366
else:
    nial = 365

obsjd = np.zeros([nfield, 10, 2])
obsmag = np.zeros([nfield, 10, 2])

obsseeing = np.zeros([nfield, 10, 2])
obsskybright = np.zeros([nfield, 10, 2])
obsairmass = np.zeros([nfield, 10, 2])

nobs1 = np.zeros(nfield)
comag1 = np.zeros(nfield)
nobs2 = np.zeros(nfield)
comag2 = np.zeros(nfield)

texp = conf['nexp']*conf['exptime'] + conf['telsettledown']

try:
    os.remove('example.log')
except:
    print('  ')
    
import time

logging.basicConfig(filename='example.log',level=logging.DEBUG)
logging.debug('This message should go to the log file')
for iday in np.arange(0, nial):
    t000 = time.time()
    obslog = []
    jd0 = tseta[iday]
    # define the filter combinations
    # ifilter = 0: ugi; 1: vrz 
    ifilter = 'vrz'
    if iday % conf['filterchangeperiod'] == 0:
        if ifilter == 'vrz':
            ifilter = 'ugi'
        else:
            ifilter = 'vrz'     
    
    # starting the sequence
    nlef = -1
    while jd0 <= trisa[iday]:             
        # achive the cloud, seeing & Sky brightness data
        days = jd2day(jd0) 
        seeing, cloud, skyzen = gw.getweather(days, tweather)

        # observable or not
        if (cloud <= conf['cloudlimit']) and (seeing <= conf['badseeinglimit']):
            vfid, vra, vdec, vam, valt, vaz, vha = selectVisible(fieldid, 
                                                     ras, decs, jd0, conf, maxAirmass=conf['maxairmass'])

            brb, distance2moon_DEG, moonAlt_DEG, brightProfile = getSkyBrightness(vra, vdec, jd0, conf, 
                                                                  extinction=0.135, 
                                                                  skyBrightness=skyzen)
            # all avaiable fields (distomoon > mindis2moon)
            indava = np.where(distance2moon_DEG>conf['mindis2moon'])
            avafid, avara, avadec = vfid[indava], vra[indava], vdec[indava]
            avaam, avaalt, avaaz, avaha = vam[indava], valt[indava], vaz[indava], vha[indava]
            avabrb = brb[indava]
            
            avafield = Table([avafid, avaam, avabrb, avaha],
                             names=['fieldid','airmass','brightness','ha'])
            
            # get the tdelay for all avafields
            if nlef ==-1:
                # the first observation
                ha0, dec0 = 0.0, conf['teldefaultdec']
                #print(avaha.size)
                tslew, tdelay = ts.teleslew(conf, ha0, dec0, nlef, avaha, avadec)
            else:
                # NOT the first observation
                obslasttime = obslog[-1]
                ra0, dec0 = obslasttime[2], obslasttime[3]
                tumy1, tumy2, ha0 = getaltaz(ra0, dec0, jd0, conf)
                tslew, tdelay = ts.teleslew(conf, ha0, dec0, nlef, avaha, avadec)                
            
            # get the avatime
            avatime = avtaa[avafid,iday]
            
            # get the obshistory
            if ifilter == 'ugi':
                avanobs, avacomag = nobs1[avafid], comag1[avafid]
            else:
                avanobs, avacomag = nobs2[avafid], comag2[avafid]
            
            # combine the avaiable information and get the ranks
            afieldinfo = gr.getfieldinfo(seeing, avafield, tdelay, avatime, avanobs, avacomag, conf)
            tranks = gr.gettotalrank(conf, afieldinfo)
            
            rfid = np.where(tranks == np.min(tranks))
            
            # get the index for the field to observe
            fid2obs = avafid[rfid[0]]
            
            # output history logs ...
            ra2obs, dec2obs = ras[fid2obs], decs[fid2obs]
            rlim2obs = (afieldinfo['rlimmag'])[rfid[0]]
            rcolim2obs = (afieldinfo['totallimitmag'])[rfid[0]]
            tdelay2obs = (afieldinfo['tdelay'])[rfid[0]]
            
            seeing2obs = seeing.copy()
            skybright2obs = (afieldinfo['brightness'])[rfid[0]]
            airmass2obs = (afieldinfo['airmass'])[rfid[0]]
            
            if ifilter == 'ugi':
                nobs1[fid2obs] += 1
                comag1[fid2obs] = rcolim2obs
                nset = nobs1[fid2obs] - 1
                nobsi = nobs1[fid2obs]
                obsjd[np.long(fid2obs[0]), np.long(nset[0]), 0] = jd0
                obsmag[np.long(fid2obs[0]), np.long(nset[0]), 0] = rlim2obs   
                obsseeing[np.long(fid2obs[0]), np.long(nset[0]), 0] = seeing2obs  
                obsskybright[np.long(fid2obs[0]), np.long(nset[0]), 0] = skybright2obs  
                obsairmass[np.long(fid2obs[0]), np.long(nset[0]), 0] = airmass2obs  

            else:
                nobs2[fid2obs] += 1
                comag2[fid2obs] = rcolim2obs            
                nset = nobs2[fid2obs] - 1
                nobsi = nobs2[fid2obs]  
                obsjd[np.long(fid2obs[0]), np.long(nset[0]), 1] = jd0
                obsmag[np.long(fid2obs[0]), np.long(nset[0]), 1] = rlim2obs                   
                obsseeing[np.long(fid2obs[0]), np.long(nset[0]), 1] = seeing2obs  
                obsskybright[np.long(fid2obs[0]), np.long(nset[0]), 1] = skybright2obs  
                obsairmass[np.long(fid2obs[0]), np.long(nset[0]), 1] = airmass2obs                 

            obslogi = [jd0, fid2obs, ra2obs, dec2obs, rlim2obs, nobsi, rcolim2obs, ifilter] 
            obslog.append(obslogi)
            logging.info('%.8f,' % obslogi[0]+ ' %i,' % obslogi[1]+ ' %.5f,' % obslogi[2]+ 
                         ' %.5f,' % obslogi[3] +' %.2f,' % obslogi[4] + ' %i,' % obslogi[5]+' %.2f.' % obslogi[6]
                        + ' ' + obslogi[7])
           
            #
            #print(t001 - t000)
            #print(tdelay2obs)
            jd0 += (tdelay2obs + texp)/3600.0/24.0
            nlef += 1
        else:
            logging.warning('Bad weather at %.8f' % jd0)
            jd0 += conf['timewait']/3600.0/24.0
            
    t001 = time.time()
    print('processing day '+ '%03i' % iday + ' taking filter combination of ' + ifilter 
          + ' with time of ' + '%.1f' % (t001-t000) + ' and fields of ' + '%i' % nlef)             

obshist = Table([obsjd, obsmag, obsseeing, obsskybright, obsairmass], 
                names=['jd','mag','seeing','skybright','airmass'])
obshist.write('obshist.fits',format='fits')
