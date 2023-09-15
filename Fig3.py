!/usr/bin/python
# -*- coding: utf-8 -*-

"""
# Author: Sentia Goursaud, University of cambridge
# Date: 06-02-2020
# Last Edited: 04-02-2020

# Description: Plotting scirpt for property-elevation relationships and lapse 
# rates (scatters) and spatial maps for high/low scenarios. Tied to 
# 'plot_spatial_100yr_clim.py', which needs to be run first to extract d18O and
# SAT data.
"""

# Execution of scripts calling module importations (eg Numpy, matplotlib,...)
execfile('/home/users/sentia/python/startup.py')

# Paths
data_root='/home/users/mholloway/my_data/bristol/ais_runs'
fig_root='/home/users/sentia/ais_paper/figures/'

# ------------------------------------------------------------------------------
# Ice core sites
# ------------------------------------------------------------------------------
# Definition (and location) of the ice core locations plot on the maps (scatter) and used in the paper
# as well as the experiment ID and names of the used simulations.

SITENAME=['VOS','DF','EDC','EDML','Taylor Dome',\
    'TALDICE','WAIS Divide','Hercules Dome','Skytrain']

lon_ice = [106.8,39.7,123.4,0.07,159.2,158.7,247.9,255,280.3]

lat_ice = [-78.5,-77.3,-75.1,-75.0,-72.8,-77.8,-79.5,-86,-78]

site_points = [('latitude', lat_ice), ('longitude', lon_ice)]
print "ice core locations:", site_points

expID=['xkgrz','xkdhi','xncmi','xncmg','xncma','xncmb','xncme','xncmf',\
    'xncmk','xncmm']

expName=['PI','BC128','BC128_1000lo','BC128_500lo','BC128_200lo',\
    'BC128_100lo','BC128_100hi','BC128_200hi','BC128_500hi',\
    'BC128_1000hi']

# ------------------------------------------------------------------------------
# load vars
# ------------------------------------------------------------------------------
# Duration for the time averaged, 50 years
tlen=50*12

# Land sea mask
filename=''+data_root+'/inidata/bbc_000cyr_mask.nc'
con = iris.Constraint(cube_func=lambda cube: cube.var_name == 'lsm')
lsm = iris.load_cube(filename, constraint=con)
lsm = lsm.collapsed(['t','surface'],iris.analysis.MEAN)
lsm.data[lsm.data == 0] = N.nan
lsm.data=N.ma.masked_invalid(lsm.data)

# d18O, SAT and orog
variables=['d18O','temp','orog']
for ID in range(len(expID)):
	name = str(expName[ID])
	print ID, ":", name
	filename=''+data_root+'/netcdf/'+expID[ID]+'_d18O_sat_z.nc'
	for VAR in variables:
		print "load vars:", str(VAR)
		con = iris.Constraint(cube_func=lambda cube: cube.var_name == VAR)
		locals()[str(name)+"_"+str(VAR)] = iris.load_cube(filename, constraint=con)
# Precip
for ID in range(len(expID)): 
    name =  str(expName[ID])
    print name
    print ID, ":", name
    root_to_data=''+data_root+'/'+expID[ID]+'/monthly'
    filename = ''+root_to_data+'/'+expID[ID]+'.precip_mm_srf.monthly.nc'
    con = iris.Constraint(cube_func=lambda cube: cube.var_name == 'precip_mm_srf')
    cube = iris.load_cube(filename, constraint=con)
    hp = cube[-tlen:].collapsed(['t','surface'],iris.analysis.MEAN)
    hp.data=N.ma.masked_invalid(hp.data)
    hp = hp*60*30*48*30.4375
    locals()[str(name)+"_P"] = hp

# ------------------------------------------------------------------------------
# orog, SAT, P, and d18O values for the ice cores 

# definition og the arrays
icecore_d18O=N.zeros([len(expID),len(lon_ice)])
icecore_temp=N.zeros([len(expID),len(lon_ice)])
icecore_P=N.zeros([len(expID),len(lon_ice)])
icecore_orog=N.zeros([len(expID),len(lon_ice)])

# nearest interporaltion for each ice core
for ID in range(len(expID)):
    cube_P = locals()[expName[ID]+"_P"]
    cube_orog = locals()[expName[ID]+"_orog"]
    cube_d18O = locals()[expName[ID]+"_d18O"]
    cube_temp = locals()[expName[ID]+"_temp"]
    for site in range(len(lon_ice)):
        points = [('latitude', lat_ice[site]), ('longitude', lon_ice[site])] 
        print(points)
        icecore_P[ID,site] = cube_P.interpolate(points, iris.analysis.Linear()).data 
        icecore_orog[ID,site] = cube_orog.interpolate(points, iris.analysis.Linear()).data 
	icecore_d18O[ID,site] = cube_d18O.interpolate(points, iris.analysis.Linear()).data 
        icecore_temp[ID,site] = cube_temp.interpolate(points, iris.analysis.Linear()).data

# Reference values: BC128, i.e. the LIG simulation
ref_orog = icecore_orog[1,:]
ref_temp=icecore_temp[1,:]
ref_d18O=icecore_d18O[1,:]
ref_P=icecore_P[1,:]

# Anomalies compared to BC128
icecore_dP=N.zeros([len(expID)-2,len(lon_ice)])
icecore_dorog=N.zeros([len(expID)-2,len(lon_ice)])
icecore_dtemp=N.zeros([len(expID)-2,len(lon_ice)])
icecore_dd18O=N.zeros([len(expID)-2,len(lon_ice)])

for ID in range(len(expID)-2):
    icecore_dP[ID]=(icecore_P[ID+2]-ref_P)/ref_P*100
    icecore_dorog[ID]=icecore_orog[ID+2]
    icecore_dd18O[ID]=(icecore_d18O[ID+2]-ref_d18O)
    icecore_dtemp[ID]=(icecore_temp[ID+2]-ref_temp)

################################################################
# ----------------------------PLOT----------------------------
################################################################

X=N.array(icecore_dorog)
Y_P=N.array(icecore_dP)
Y_temp=N.array(icecore_dtemp)
Y_dd18O =N.array(icecore_dd18O)

filname='AIS_FIG2_dXvsdZ_1km_refBC128_30yr_Feb20' # Filename to save the figure 
colors=['c','RoyalBlue','Indigo','m','Grey','k','Orange','Crimson','green'] # colors associated at each ice core
markertype=['o','s','o','s','o','s','o','s','o']
site_name=SITENAME
Nsite = len(site_name)

fig = plt.figure(figsize=[8,10])
#---------------------------
ax1 = fig.add_subplot(3,1,2)
ax1.set_ylabel(r'$\frac{\Delta{P}}{P_{Ref}}$ (%)',fontsize=14)    
for SITE in range(Nsite):
    ax1.plot([X.min()-500,X.max()+500],[0,0],':k')
    ax1.scatter(X[:,SITE],Y_P[:,SITE], edgecolor=colors[SITE], color='White',marker=markertype[SITE],s=50)
# polynomial fit
for i in range(Nsite-1):
    a,b,c = np.polyfit(X[:,i],Y_P[:,i],2)
    y_poly2 = c+b*X[:,i]+a*X[:,i]**2
    ax1.plot(X[:,i],y_poly2,color=colors[i],linestyle='--',lw=1.0,alpha=0.5)
# Skytrain linear fit
a,b = np.polyfit(X[:,-1],Y_P[:,-1],1)
y_sky = b+a*X[:,-1]
ax1.plot(X[:,-1],y_sky,color='g',linestyle='-',lw=1.0,alpha=0.5)
gca().set_xticklabels(['']*6)
ax1.annotate('B', xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
#---------------------------
ax2 = fig.add_subplot(3,1,1)
ax2.set_ylabel(r'$\Delta{SAT}$ ($^{\circ}$C)',fontsize=12)
for SITE in range(Nsite):
    ax2.plot([X.min()-500,X.max()+500],[0,0],':k')
    ax2.scatter(X[:,SITE],Y_temp[:,SITE], edgecolor=colors[SITE], color='White',marker=markertype[SITE],s=50)
# linear fit
for i in range(Nsite-1):
    a,b = np.polyfit(X[:,i],Y_temp[:,i],1)
    y_lin = b+a*X[:,i]
    ax2.plot(X[:,i],y_lin,color=colors[i],linestyle='-',lw=1.0,alpha=0.5)
# Skytrain fit
a,b,c = np.polyfit(X[:,-1],Y_temp[:,-1],2)
y_sky = c+b*X[:,-1]+a*X[:,-1]**2
ax2.plot(X[:,-1],y_sky,color='g',linestyle='--',lw=1.0,alpha=0.5)
gca().set_xticklabels(['']*6)
ax2.annotate('A', xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
#---------------------------
ax3 = fig.add_subplot(3,1,3)
ax3.set_xlabel('Ice core site elevation (m)',fontsize=12)
ax3.set_ylabel(r'$\Delta\delta^{18}$O (permille)',fontsize=12)   
#ax3.set_ylabel(r'$\Delta{P}$ (mm month$^{-1}$)',fontsize=12)#^{18}$O (permille)',fontsize=12)    
for SITE in range(Nsite):
    ax3.plot([X.min()-500,X.max()+500],[0,0],':k')
    ax3.scatter(X[:,SITE],Y_dd18O[:,SITE], edgecolor=colors[SITE], color='White',marker=markertype[SITE],s=50)
# regression fits
for i in np.arange(5,9):
    a,b = np.polyfit(X[:,i],Y_dd18O[:,i],1)
    y_lin = b+a*X[:,i]
    ax3.plot(X[:,i],y_lin,color=colors[i],linestyle='-',lw=1.0,alpha=0.5)
# Polynomial fits
for i in np.arange(Nsite-4):
    a,b,c = np.polyfit(X[:,i],Y_dd18O[:,i],2)
    y_poly2 = c+b*X[:,i]+a*X[:,i]**2
    ax3.plot(X[:,i],y_poly2,color=colors[i],linestyle='--',lw=1.0,alpha=0.5)
ax3.annotate('C', xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')

xtext=4100
ytext=11
dytxt=0.5
dxtxt=0.05
for kk in range(Nsite):
    ax2.scatter(xtext, ytext-dytxt,c='White',s=30, marker=markertype[kk], edgecolors=colors[kk])
    ax2.text(xtext+100, ytext-dytxt-(dxtxt), str(site_name[kk]),fontsize=11,color='k')
    dytxt=dytxt+2

plt.tight_layout()
plt.savefig(''+fig_root+'/'+filname+'.png', dpi=300)
plt.savefig(''+fig_root+'/'+filname+'.pdf', dpi=300)

#END

