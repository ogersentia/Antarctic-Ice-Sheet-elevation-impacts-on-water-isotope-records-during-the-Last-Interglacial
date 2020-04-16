!/usr/bin/python
# -*- coding: utf-8 -*-

"""
# Author: Sentia Goursaud, University of cambridge
# Date: 02-03-2020
# Last Edited: 02-03-2020

# Description: Plotting maps of elevational relationships.
"""

# Exceution of scripts calling module importations (eg Numpy, matplotlib,...)
execfile('/home/users/mholloway/python/startup.py')

# Paths 
data_root='/home/users/mholloway/my_data/bristol/ais_runs'
fig_root='/home/users/sentia/ais_paper/figures/'

# ------------------------------------------------------------------------------
# Ice core sites
# ------------------------------------------------------------------------------
# Definition (and location) of the ice core locations plot on the maps (scatter) and used in the paper,
# as well as the experiment ID and names of the used simulations.
SITENAME=['VOS','DF','EDC','EDML','TALDICE','Taylor Dome','WAIS Divide','Hercules Dome','Skytrain']

lon_ice = [106.8,39.7,123.4,0.07,159.2,158.7,247.9,255,280.3]

lat_ice = [-78.5,-77.3,-75.1,-75.0,-72.8,-77.8,-79.5,-86,-78]

site_points = [('latitude', lat_ice), ('longitude', lon_ice)]
print "ice core locations:", site_points

expID=['xkgrz','xkdhi','xncmi', 'xncmg','xncma','xncmb','xncme','xncmf','xncmk','xncmm']

expName=['PI','BC128','BC128_1000lo','BC128_500lo','BC128_200lo','BC128_100lo','BC128_100hi',\
        'BC128_200hi','BC128_500hi','BC128_1000hi']

# ------------------------------------------------------------------------------
# load variables
# ------------------------------------------------------------------------------
tlen=50*12 #Duration of the time-averaged

# Land sea mark
filename=''+data_root+'/inidata/bbc_000cyr_mask.nc'
con = iris.Constraint(cube_func=lambda cube: cube.var_name == 'lsm')
LSM = iris.load_cube(filename, constraint=con)
lsm = LSM.collapsed(['t','surface'],iris.analysis.MEAN)
lsm.data[lsm.data == 0] = np.nan
lsm.data=np.ma.masked_invalid(lsm.data)

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

# Arrays of simulations used in the analysis
IDs=['PI','BC128_1000lo','BC128_500lo','BC128_200lo','BC128_100lo',
     'BC128_500hi','BC128_1000hi','BC128_200hi','BC128_100hi']
orogvars=[]
satvars=[]
Pvars=[]
isovars=[]
for ii in [expName[expName.index(x)] for x in IDs]:
    orogvars.append(locals()[str(ii)+"_orog"])
    satvars.append(locals()[str(ii)+"_temp"])
    Pvars.append(locals()[str(ii)+"_P"])
    isovars.append(locals()[str(ii)+"_d18O"])    

# References _ Present days
reforog=PI_orog
reftemp=PI_temp
refprecip=PI_P
refd18O=PI_d18O

# ------------------------------------------------------------------------------
# Spatial linear correlations -- FOR MAPS
# ------------------------------------------------------------------------------
land_ix = np.argwhere(lsm.data==1)

p_temp=np.zeros((73,96))
a_temp=np.zeros((73,96))
r_temp=np.zeros((73,96))
for i in range(len(land_ix)):
    lat_ix = land_ix[i,0]
    lon_ix = land_ix[i,1]
    x=[orogvars[ID].data[lat_ix,lon_ix] -reforog.data[lat_ix,lon_ix] for ID in range(len(IDs))]
    y=[satvars[ID].data[lat_ix,lon_ix] -reftemp.data[lat_ix,lon_ix] for ID in range(len(IDs))]
    a_temp[lat_ix,lon_ix],b_temp,r_temp[lat_ix,lon_ix],p_temp[lat_ix,lon_ix],stderr=stats.mstats.linregress(x,y)

p_precip=np.zeros((73,96))
a_precip=np.zeros((73,96))
r_precip=np.zeros((73,96))
for i in range(len(land_ix)):
    lat_ix = land_ix[i,0]
    lon_ix = land_ix[i,1]
    x=[orogvars[ID].data[lat_ix,lon_ix] -reforog.data[lat_ix,lon_ix] for ID in range(len(IDs))]
    y=[Pvars[ID].data[lat_ix,lon_ix] -refprecip.data[lat_ix,lon_ix] for ID in range(len(IDs))]
    a_precip[lat_ix,lon_ix],b_precip,r_precip[lat_ix,lon_ix],p_precip[lat_ix,lon_ix],stderr=stats.mstats.linregress(x,y)

p_d18O=np.zeros((73,96))
a_d18O=np.zeros((73,96))
r_d18O=np.zeros((73,96))
for i in range(len(land_ix)):
    lat_ix = land_ix[i,0]
    lon_ix = land_ix[i,1]
    x=[orogvars[ID].data[lat_ix,lon_ix] -reforog.data[lat_ix,lon_ix] for ID in range(len(IDs))]
    y=[isovars[ID].data[lat_ix,lon_ix] -refd18O.data[lat_ix,lon_ix] for ID in range(len(IDs))]
    a_d18O[lat_ix,lon_ix],b_d18O,r_d18O[lat_ix,lon_ix],p_d18O[lat_ix,lon_ix],stderr=stats.mstats.linregress(x,y)

# Statistical outputs multiplied by the land sea mask, and slope in /100m.
a_temp = a_temp*lsm.data*100
r_temp = r_temp*lsm.data
p_temp = p_temp*lsm.data
a_precip = a_precip*lsm.data*100
r_precip = r_precip*lsm.data
p_precip = p_precip*lsm.data
a_d18O = a_d18O*lsm.data*100
r_d18O = r_d18O*lsm.data
p_d18O = p_d18O*lsm.data

# ---------------------------------------------------------------
# PLOT 
# ---------------------------------------------------------------
# colormaps
import string
import matplotlib.colors as mcolors
colors1 = plt.cm.YlOrBr_r(np.linspace(0., 1, 128))
colors2 = plt.cm.BuPu(np.linspace(0, 1, 128))
colors = np.vstack((colors1, colors2))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

red = N.array([0, 0, 0, 10, 0, 60, 72, 221, 255,255, 250, 255, 229, 239, 228, 205, 161, 116, 77]) / 256.
green = N.array([0, 0, 66, 144, 191, 209, 255, 242, 255,255, 250, 255, 235, 190, 128, 72, 33, 29, 30]) / 256.
blue = N.array([128, 255, 255, 255, 255, 204, 250, 243, 255, 255, 210, 0, 99, 63, 39, 27, 22, 29, 27]) / 256.
colormap = N.array([red, green, blue]).T
my_colormap = matplotlib.colors.ListedColormap(colormap, name='my_colormap')
fmt='%.0f'

# latitude and longitude
lat,lon=reftemp.coord('latitude').points,reftemp.coord('longitude').points
lon1=reftemp.coord('longitude').points

fig = plt.figure(figsize=[8,8])
#---------------------------
ax1 = fig.add_subplot(2,3,1)
map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
map.drawcoastlines()
map.drawcountries()
map.drawparallels(N.arange(-80.,81.,20.), color='grey')
map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
var,lonsout = mpl_toolkits.basemap.shiftgrid(180,a_temp,lon,start=False, cyclic=360.0)
var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
p,lonsout1 = mpl_toolkits.basemap.shiftgrid(180,p_temp,lon1,start=False, cyclic=360.0)
p,lonsout1 = mpl_toolkits.basemap.addcyclic(p,lonsout1)
x,y = map(*np.meshgrid(lonsout,lat))
cs = map.contourf(x,y,var,cmap=plt.cm.rainbow,alpha=0.7,levels=np.arange(-1.5,1.6,0.1))
#cbar.set_label(r'${\Delta}$Z-${\Delta}$SAT slope (${\circ}$C / 100m )',fontsize=12)#${\Delta}$Z-${\Delta}$SAT
xpt,ypt = map(lon_ice,lat_ice)
lonpt, latpt = map(xpt,ypt,inverse=True)
collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
cs2=map.contourf(x,y,p,levels=[0,0.05],colors='none',hatches=['//',''],extend='both')
cs3=map.contour(x,y,p,colors='k',levels=[0.05])
ax1.annotate('A', xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
ax1.set_title(r'$\Delta$SAT (${\circ}$C)')
colorbar_axes = fig.add_axes([0.08, 0.6, 0.005, 0.2])
cbar = plt.colorbar(cs, orientation= 'vertical',cax=colorbar_axes)
cbar.set_ticks(N.arange(-1.5,1.6,0.5))
cbar.set_ticklabels([N.round(x,1) for x in N.arange(-1.5,1.6,0.5)])
cbar.set_label(r'Slope (/100m)',fontsize=12)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().yaxis.set_label_position('left')

#---------------------------#
ax2 =  fig.add_subplot(2,3,4)
map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
map.drawcoastlines()
map.drawcountries()
map.drawparallels(N.arange(-80.,81.,20.), color='grey')
map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
var,lonsout = mpl_toolkits.basemap.shiftgrid(180,r_temp,lon,start=False, cyclic=360.0)
var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
x,y = map(*np.meshgrid(lonsout,lat))
cs = map.contourf(x,y,var*var,cmap=plt.cm.rainbow,alpha=0.7,levels=np.arange(0,1.1,0.1))
#cbar = plt.colorbar(cs, orientation= 'horizontal')
#cbar.set_label(r'${\Delta}$Z-${\Delta}$SAT correlation coefficient',fontsize=12)#${\Delta}$Z-${\Delta}$SAT
collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
cs2=map.contourf(x,y,p,levels=[0,0.05],colors='none',hatches=['//',''],extend='both')
cs3=map.contour(x,y,p,colors='k',levels=[0.05])
ax2.annotate('B', xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
colorbar_axes = fig.add_axes([0.08, 0.2, 0.005, 0.2])
cbar = plt.colorbar(cs, orientation= 'vertical',cax=colorbar_axes)
cbar.set_ticks(N.arange(0,1.2,0.2))
cbar.set_ticklabels([N.round(x,1) for x in N.arange(0,1.2,0.2)])
cbar.set_label(r'r$^{2}$',fontsize=12) #Correlation coefficient
plt.gca().yaxis.set_ticks_position('left')
plt.gca().yaxis.set_label_position('left')

#---------------------------
ax3 = fig.add_subplot(2,3,2)
map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
map.drawcoastlines()
map.drawcountries()
map.drawparallels(N.arange(-80.,81.,20.), color='grey')
map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
var,lonsout = mpl_toolkits.basemap.shiftgrid(180,a_precip,lon,start=False, cyclic=360.0)
var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
p,lonsout1 = mpl_toolkits.basemap.shiftgrid(180,p_precip,lon1,start=False, cyclic=360.0)
p,lonsout1 = mpl_toolkits.basemap.addcyclic(p,lonsout1)
x,y = map(*np.meshgrid(lonsout,lat))
cs = map.contourf(x,y,var,cmap=plt.cm.rainbow,alpha=0.7,levels=np.arange(-1.5,1.6,0.1))
#cbar = plt.colorbar(cs, orientation= 'vertical')
#cbar.set_label(r'${\Delta}$Z-${\Delta}$P slope (mm?month / 100m )',fontsize=12)#${\Delta}$Z-${\Delta}$SAT
xpt,ypt = map(lon_ice,lat_ice)
lonpt, latpt = map(xpt,ypt,inverse=True)
collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
cs2=map.contourf(x,y,p,levels=[0,0.05],colors='none',hatches=['//',''],extend='both')
cs3=map.contour(x,y,p,colors='k',levels=[0.05])
cs4=map.contour(x,y,var,colors='darkblue',levels=[-50,-20,-5])
ax3.clabel(cs4, cs4.levels,inline=True, fontsize=10)
ax3.annotate('C', xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
ax3.set_title('$\Delta$P (mm/month)')

#---------------------------
ax4 =  fig.add_subplot(2,3,5)
map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
map.drawcoastlines()
map.drawcountries()
map.drawparallels(N.arange(-80.,81.,20.), color='grey')
map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
var,lonsout = mpl_toolkits.basemap.shiftgrid(180,r_precip,lon,start=False, cyclic=360.0)
var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
cs = map.contourf(x,y,var*var,cmap=plt.cm.rainbow,alpha=0.7,levels=np.arange(0,1.1,0.1))
#cbar = plt.colorbar(cs, orientation= 'horizontal')
#cbar.set_label(r'${\Delta}$Z-${\Delta}$P correlation coefficient',fontsize=12)#${\Delta}$Z-${\Delta}$SAT
collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
cs2=map.contourf(x,y,p,levels=[0,0.05],colors='none',hatches=['//',''],extend='both')
cs3=map.contour(x,y,p,colors='k',levels=[0.05])
ax4.annotate('D', xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')

#---------------------------
ax5 = fig.add_subplot(2,3,3)
map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
map.drawcoastlines()
map.drawcountries()
map.drawparallels(N.arange(-80.,81.,20.), color='grey')
map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
var,lonsout = mpl_toolkits.basemap.shiftgrid(180,a_d18O,lon,start=False, cyclic=360.0)
var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
p,lonsout1 = mpl_toolkits.basemap.shiftgrid(180,p_d18O,lon1,start=False, cyclic=360.0)
p,lonsout1 = mpl_toolkits.basemap.addcyclic(p,lonsout1)
x,y = map(*np.meshgrid(lonsout,lat))
cs = map.contourf(x,y,var,cmap=plt.cm.rainbow,alpha=0.7,levels=np.arange(-1.5,1.6,0.1))
#cbar = plt.colorbar(cs, orientation= 'horizontal')
#cbar.set_label(r'${\Delta}$Z-$\Delta\delta^{18}$O slope ('+u' (\u2030)'+'/ 100m )',fontsize=12)#${\Delta}$Z-${\Delta}$SAT
xpt,ypt = map(lon_ice,lat_ice)
lonpt, latpt = map(xpt,ypt,inverse=True)
collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
cs2=map.contourf(x,y,p,levels=[0,0.05],colors='none',hatches=['//',''],extend='both')
cs3=map.contour(x,y,p,colors='k',levels=[0.05])
cs4=map.contour(x,y,var,colors='darkblue',levels=[-50,-20,-5])
ax5.clabel(cs4, cs4.levels,inline=True, fontsize=10)
ax5.annotate('E', xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
ax5.set_title('$\Delta\delta^{18}$O '+u'(\u2030)' )
#---------------------------
ax6 =  fig.add_subplot(2,3,6)
map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
map.drawcoastlines()
map.drawcountries()
map.drawparallels(N.arange(-80.,81.,20.), color='grey')
map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
var,lonsout = mpl_toolkits.basemap.shiftgrid(180,r_d18O,lon,start=False, cyclic=360.0)
var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
cs = map.contourf(x,y,var*var,cmap=plt.cm.rainbow,alpha=0.7)#,levels=np.arange(0,1.1,0.1))
#cbar = plt.colorbar(cs, orientation= 'horizontal')
#cbar.set_label(r'${\Delta}$Z-$\Delta\delta^{18}$O correlation coefficient',fontsize=12)#${\Delta}$Z-${\Delta}$SAT
collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
cs2=map.contourf(x,y,p,levels=[0,0.05],colors='none',hatches=['//',''],extend='both')
cs3=map.contour(x,y,p,colors='k',levels=[0.05])
ax6.annotate('F', xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')

#plt.tight_layout()

# ----------------------------------------------------------
filname='map_elevation_relationships'
plt.savefig(''+fig_root+filname+'.png', dpi=300)
plt.savefig(''+fig_root+filname+'.pdf', dpi=300)

#END
