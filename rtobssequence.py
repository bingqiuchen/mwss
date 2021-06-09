import matplotlib.pyplot as plt 
from astropy.table import Table, Column
from astropy.io import fits
import getrealweather as gw
import getrank as gr
import telescope as ts
import dithering as dt
from sitesky import *
import numpy as np
import time as tm
from tqdm import trange
import logging

"""
revised in Mar 31 2021, simulate the real case; using the real time weather
          obtained from the Lijiang Obnservatory
revised in Jan 12 2021, adding the dithering processes
          0.06<d𝑅𝐴=d𝐷𝑒𝑐<0.1 deg
revised in Aug 31 2020, adding additional outputs
          such as airmass, v-band sky brightness and seeing
Written by BQC in 201912
"""

def timewaitshow(tsleep):
    """
    Show the waiting time tsleep remaining with bars
    """
    Nshow = int(np.ceil(tsleep))
    for i in trange(Nshow):
        tm.sleep(1)
        
def timetoshow(tarjd):
    """
    Show the time to targeting jd with remaining of bars
    """
    nowjd = gw.getjdn()
    tsleep = (tarjd - nowjd) * 24.0 * 3600.0
    Nshow = int(np.ceil(tsleep))
    logging.info('%.1f' %tsleep + ' seconds left to JD %.8f' % tarjd)
    for i in trange(Nshow):
        tm.sleep(1)

try:
    os.remove('example.log')
except:
    print('  ')
    
logging.basicConfig(filename='example.log',level=logging.DEBUG)
logging.debug('This message should go to the log file')

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

if isLeap(conf['year']):
    nial = 366
else:
    nial = 365

try: 
    with open('test.npy', 'rb') as f:
        obsjd       =np.load(f)
        obsmag      =np.load(f)
        obsdith     =np.load(f)
        obsseeing   =np.load(f)
        obsskybright=np.load(f)
        obsairmass  =np.load(f)
        nobs1       =np.load(f)
        comag1      =np.load(f)
        nobs2       =np.load(f)
        comag2      =np.load(f)
except:
    obsjd = np.zeros([nfield, 10, 2])
    obsmag = np.zeros([nfield, 10, 2])
    obsdith = np.zeros([nfield, 10, 2])
   
    obsseeing = np.zeros([nfield, 10, 2])
    obsskybright = np.zeros([nfield, 10, 2])
    obsairmass = np.zeros([nfield, 10, 2])
   
    nobs1 = np.zeros(nfield)
    comag1 = np.zeros(nfield)
    nobs2 = np.zeros(nfield)
    comag2 = np.zeros(nfield)

texp = conf['nexp']*conf['exptime'] + conf['telsettledown']

# run real time observation plan selection
jdn = gw.getjdn()
nlef = -1
obslog = []
deltaradec = np.random.uniform(conf['mindithering'],conf['maxdithering'],1500)  # dithering
while (jdn < 2459518):
    # get present time and if in night then read weather data
    dayn = gw.getdayn()
    jdn = gw.getjdn()

    if (jdn <= tseta[dayn]) | (jdn >= trisa[dayn]):
    #if (jdn >= tseta[dayn]) & (jdn <= trisa[dayn]):
        nlef = -1
        obslog = []
        deltaradec = np.random.uniform(conf['mindithering'],conf['maxdithering'],1500)  # dithering
        
        # determine the sun set time
        tsetime = tseta[dayn]
        if tsetime < jdn:
          tsetime = tseta[dayn+1]

        tsleep = (tsetime - jdn) * 24.0 * 3600.0
        logging.warning('It is daytime now, MSSS will sleep %i seconds untill Sunset!' % tsleep)
        #timewaitshow(tsleep)
        timetoshow(tsetime)
        logging.info('Night comes, MSSS now start for observations.')
    else:
        tc, cloud, tsee, seeing, td, skybright = gw.getweather()
        #for test we will create time manually
        #timenow = tm.strftime("%Y-%m-%dT%H:%M:%S",tm.gmtime()) 
        #tnow = Time(str(timenow), format='isot', scale='utc')
        #tc, cloud, tsee, seeing, td, skybright = tnow, 0, tnow, 1.28, tnow, 21.5

        # check if the weather data collected recently 
        if (np.abs((jdn-tc.jd)*24.0*60.0) > 10.0) | (np.abs((jdn-td.jd)*24.0*60.0) > 10.0) | (np.abs((jdn-tsee.jd)*24.0*60.0) > 10.0):
        #if (np.abs((jdn-tc.jd)*24.0*60.0) > 5.0) | (np.abs((jdn-td.jd)*24.0*60.0) > 5.0):
            #breakpoint()
            logging.warning('''Weather data too old, please check! 
                Cloud time: ''' +tc.strftime("%Y/%m/%dT%H:%M:%S")+'''; 
                Seeing time: '''+tsee.strftime("%Y/%m/%dT%H:%M:%S")+'''; 
                Sky time: '''+td.strftime("%Y/%m/%dT%H:%M:%S"))
            twaitw = conf['timewait']
            timewaitshow(twaitw)
        else:
            # check the weather data
            #skybright = np.array(21.0)
            if (cloud > conf['cloudlimit']) | (seeing > conf['badseeinglimit']) | (skybright < conf['skybrightlimit']):
                logging.warning('Bad weather at JD %.8f: ' % jdn + 
                      '; with seeing: %.2f ' % seeing + ', skybrightness: %.1f' % 
                      skybright + ' and cloud: %i.' % cloud)
                twaitw = conf['timewait']
                logging.warning('waiting for good weather, time remaining: %i seconds.' % twaitw)
                timewaitshow(twaitw)
            else:
                # cloud, seeing, skybright = np.array(cloud), np.array(seeing), np.array(skybright)
                # define the filter combinations
                # ifilter = 0: ugi; 1: vrz 
                ifilter = 'vrz'
                if dayn % conf['filterchangeperiod'] == 0:
                    if ifilter == 'vrz':
                        ifilter = 'ugi'
                        #logging.info('using the filter combinations of ugi in day %i' % dayn)
                    else:
                        ifilter = 'vrz'   
                        #logging.info('using the filter combinations of vrz in day %i' % dayn)

                # make observations
                vfid, vra, vdec, vam, valt, vaz, vha = selectVisible(fieldid, 
                                                         ras, decs, jdn, conf, maxAirmass=conf['maxairmass'])

                brb, distance2moon_DEG, moonAlt_DEG, brightProfile = getSkyBrightness(vra, vdec, jdn, conf, 
                                                                      extinction=0.135, 
                                                                      skyBrightness=skybright)
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
                    #logging.info(avaha.size)
                    tslew, tdelay = ts.teleslew(conf, ha0, dec0, nlef, avaha, avadec)
                else:
                    # NOT the first observation
                    obslasttime = obslog[-1]
                    ra0, dec0 = obslasttime[2], obslasttime[3]
                    tumy1, tumy2, ha0 = getaltaz(ra0, dec0, jdn, conf)
                    tslew, tdelay = ts.teleslew(conf, ha0, dec0, nlef, avaha, avadec)                

                # get the avatime
                avatime = avtaa[avafid,dayn]

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

                # dithering
                dtradec = deltaradec[nlef+1]

                if ifilter == 'ugi':
                    ra2obs1, dec2obs1, ra2obs2, dec2obs2 = dt.dithering(ra2obs,dec2obs,nobs1[fid2obs],dtradec)
                    nobs1[fid2obs] += 1
                    comag1[fid2obs] = rcolim2obs
                    nset = nobs1[fid2obs] - 1
                    nobsi = nobs1[fid2obs]
                    obsjd[np.long(fid2obs[0]), np.long(nset[0]), 0] = jdn
                    obsmag[np.long(fid2obs[0]), np.long(nset[0]), 0] = rlim2obs   
                    obsdith[np.long(fid2obs[0]), np.long(nset[0]), 0] = dtradec 
                    obsseeing[np.long(fid2obs[0]), np.long(nset[0]), 0] = seeing2obs  
                    obsskybright[np.long(fid2obs[0]), np.long(nset[0]), 0] = skybright2obs  
                    obsairmass[np.long(fid2obs[0]), np.long(nset[0]), 0] = airmass2obs  

                else:
                    ra2obs1, dec2obs1, ra2obs2, dec2obs2 = dt.dithering(ra2obs,dec2obs,nobs2[fid2obs],dtradec)
                    nobs2[fid2obs] += 1
                    comag2[fid2obs] = rcolim2obs            
                    nset = nobs2[fid2obs] - 1
                    nobsi = nobs2[fid2obs]  
                    obsjd[np.long(fid2obs[0]), np.long(nset[0]), 1] = jdn
                    obsmag[np.long(fid2obs[0]), np.long(nset[0]), 1] = rlim2obs  
                    obsdith[np.long(fid2obs[0]), np.long(nset[0]), 1] = dtradec  
                    obsseeing[np.long(fid2obs[0]), np.long(nset[0]), 1] = seeing2obs  
                    obsskybright[np.long(fid2obs[0]), np.long(nset[0]), 1] = skybright2obs  
                    obsairmass[np.long(fid2obs[0]), np.long(nset[0]), 1] = airmass2obs                 

                obslogi = [jdn, fid2obs, ra2obs1, dec2obs1, ra2obs2, dec2obs2, rlim2obs, nobsi, rcolim2obs, ifilter] 
                obslog.append(obslogi)
                logging.info('Observation plan: %.8f,' % obslogi[0]+ ' %i,' % obslogi[1]+ ' %.5f,' % obslogi[2]+ 
                             ' %.5f,' % obslogi[3] + ' %.5f,' % obslogi[4]+ 
                             ' %.5f,' % obslogi[5] +' %.2f,' % obslogi[6] + ' %i,' % obslogi[7]+' %.2f.' % obslogi[8]
                            + ' ' + obslogi[9])
               
                nlef += 1

                # save files:
                with open('test.npy', 'wb') as f:
                    np.save(f,obsjd)
                    np.save(f,obsmag)
                    np.save(f,obsdith)
                    np.save(f,obsseeing)
                    np.save(f,obsskybright)
                    np.save(f,obsairmass)
                    np.save(f,nobs1)
                    np.save(f,comag1)
                    np.save(f,nobs2)
                    np.save(f,comag2)

                #jdnow = gw.getjdn()
                #twait = (tdelay2obs + texp) - ((jdnow-jdn)*3600.0)
                logging.info('Waiting the telescope to finish the current visit.')
                jdtar = jdn + ((tdelay2obs + texp)/3600.0/24.0)
                timetoshow(jdtar)
                #timewaitshow(twait)
                 

