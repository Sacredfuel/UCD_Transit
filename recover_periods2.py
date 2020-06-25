import everest
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.timeseries import BoxLeastSquares
import pandas as pd
import numpy as np
import lightkurve as lk
import math
plt.rc('text', usetex=False)
sns.set_context('talk')

df = pd.DataFrame(columns = ['EPIC','PERIOD','DEPTH','TRANSITTIME','Rp/Rs'])
# how to 'SpT','Rs','Rp'
#['a', 'b', 'c', 'd', 'e']
#'EPIC','PERIOD','DEPTH','TRANSITTIME','Rp/Rs'

def currentFrame():
    return df

class model:
    def __init__(self, EPIC, K2name, minperiod = 1, maxperiod = 30, window_length=301, sigma=10, addToFrame = True):
        self.K2name = K2name
        self.EPIC = EPIC
        for i in np.arange(1,18,1):
            try:
                self.star = everest.Everest(EPIC, quiet=True, season=i)
            except Exception as e:
                i+=i

        self.rawlc = lk.LightCurve(time = self.star.time, flux = self.star.flux)

        self.flat = self.rawlc.flatten(window_length=window_length).remove_outliers(sigma=sigma)

        model = BoxLeastSquares(self.flat.time, self.flat.flux, dy=0.01)
        testperiods = np.arange(minperiod, maxperiod, 0.001)
        self.periodogram = model.power(testperiods, 0.16)
        maxID = np.argmax(self.periodogram.power)
        self.best_fit_period = self.periodogram.period[maxID]
        self.transit_time = self.periodogram.transit_time[maxID]
        self.transit_depth = self.periodogram.depth[maxID]

        self.foldedlc = self.flat.fold(period = self.best_fit_period, t0 = self.transit_time)

        if(addToFrame==True):
            df.append([self.EPIC, self.best_fit_period, self.transit_depth, self.transit_time, math.sqrt(self.transit_depth)])


    def dvs(self):
        return self.star.dvs()
    def plotRaw(self):
        fig = plt.figure(figsize=(10,4))
        plt.plot(self.rawlc.time, self.rawlc.flux,marker='o',markersize=2, linestyle='None',color='k')
        plt.xlabel('MJD')
        plt.ylabel('Flux')
        plt.savefig('rawLC_'+self.K2name+'.jpg',bbox_extra_artists=(),bbox_inches='tight')
    def plotFlat(self):
        fig = plt.figure(figsize=(10,4))
        plt.plot(self.flat.time, self.flat.flux,marker='o',markersize=2, linestyle='None',color='k')
        #plt.ylim(percentile(self.flat.flux, 5), percentile(self.flat.flux, 95))
        plt.xlabel('MJD')
        plt.ylabel('Flux')
        plt.savefig('flatLC_'+self.K2name+'.jpg',bbox_extra_artists=(),bbox_inches='tight')
    def getPeriods(self):
        return self.best_fit_period, self.transit_time, self.transit_depth
    def periodogramPlot(self):
        plt.plot(self.periodogram.period, self.periodogram.power)
        plt.xlabel('Period (days)')
        plt.ylabel('Power')
        plt.savefig('periodogram_'+self.K2name+'.jpg')
        return periodogramPlot
    def plotFolded(self):
        fig = plt.figure(figsize=(10,4))
        plt.plot(self.foldedlc.time, self.foldedlc.flux,marker='o',markersize=2,linestyle='None',color='k')
        plt.xlabel('Phase')
        plt.ylabel('Flux')
        plt.savefig('foldedLC_'+self.K2name+'.jpg',bbox_extra_artists=(),bbox_inches='tight')
    def getRadii(self):
        Rsun = 6.955 * 10**10 #cm
        REarth = 6.378*10**8 #cm
        df = pd.DataFrame(index = np.arange(12), columns = ('EPIC','PERIOD','DEPTH','Rp/Rs','SpT','Rs','Rp'))
        df['EPIC'] = self.EPIC
        df['SpT']  = ['M1','M3','M2'] 
        rsdict = {'M0':0.62, 'M1':0.49, 'M2':0.44, 'M3':0.39, 'N/A':0}
        df['PERIOD'] = self.best_fit_period
        df['DEPTH'] = self.transit_depth
        df['Rs'] = df['SpT'].map(rsdict)*Rsun #cm
        df['Rp/Rs'] = np.sqrt(transit_depth)
        df['Rp'] = df['Rs']*df['Rp/Rs'] / REarth
        return df

