!/usr/bin/python
# -*- coding: utf-8 -*-


"""
# Author: Max Holloway, BAS 
# Date: 02-2018
# Modifications by Sentia Goursaud, University of Cambridge
# Last Edited: 09-03-2020

# Description: Plotting Spatial maps for Antarctic Ice Sheet high/low scenarios run by the HadCM3 model with 128 ka orbital parameters - for the orography with the sea ice extent, the precipitation, SAT and  d18O with the mean sea level pressures.  
"""

# Exceution of scripts calling module importations (eg Numpy, matplotlib,...)
execfile('/home/users/sentia/python/startup.py')

# Paths 
data_root='/home/users/mholloway/my_data/bristol/ais_runs' #Root path for the original simulations
fig_root='/home/users/sentia/ais_paper/figures/' #For saving the figures


# ------------------------------------------------------------------------------
# Ice core sites and experiments def.
# ------------------------------------------------------------------------------
# Definition (and location) of the ice core locations plot on the maps (scatter) and used in the paper,
# as well as the experiment ID and names of the used simulations.

SITENAME=['VOS','DF','EDC','EDML','TALDICE','Taylor Dome','WAIS Divide',\
    'Hercules Dome','Skytrain']

lon_ice = [106.8,39.7,123.4,0.07,159.2,158.7,247.9,255,280.3]

lat_ice = [-78.5,-77.3,-75.1,-75.0,-72.8,-77.8,-79.5,-86,-78]

expID=['xkgrz','xkdhi','xncmi','xncmg','xncma','xncmb','xncme','xncmf',\
    'xncmk','xncmm']

expName=['PI','BC128','BC128_1000lo','BC128_500lo','BC128_200lo',\
    'BC128_100lo','BC128_100hi','BC128_200hi','BC128_500hi',\
    'BC128_1000hi']

# ---------------------------------------------------------------------------------------
# LOAD VARS
# ---------------------------------------------------------------------------------------
#Load of the variables using iris (https://scitools.org.uk/iris/docs/latest/index.html)

# 1. Land sea mask
filename=''+data_root_max+'/inidata/bbc_000cyr_mask.nc'
con = iris.Constraint(cube_func=lambda cube: cube.var_name == 'lsm')
LSM = iris.load_cube(filename, constraint=con)
lsm = LSM.collapsed(['t','surface'],iris.analysis.MEAN)
lsm.data[lsm.data == 0] = np.nan
lsm.data=np.ma.masked_invalid(lsm.data)

# 2. sea ice
# Function to exctract September sea-ice
def extract_sep(var):
    iris.coord_categorisation.add_month(var, 't', name='month')
    month = var.coord('month')
    sep = iris.Constraint(month=['Sep'])
    cube_sep = var.extract(sep)
    print repr(cube_sep)
    return month, cube_sep

tlen=50*12 # duration for the time average, 50yrs
for ID in range(len(expID)):
    name = str(expName[ID])
    print ID, ":", name
    filename = ''+data_root_max+'/'+expID[ID]+'/monthly/'+expID[ID]+'.iceconc_mm_srf.monthly.nc'
    con = iris.Constraint(cube_func=lambda cube: cube.var_name == 'iceconc_mm_srf')
    cube  = iris.load_cube(filename, constraint=con)
    ice = cube[-tlen:].collapsed('surface',iris.analysis.MEAN)
    month, ice_sep = extract_sep(ice)
    locals()[str(name)+"_sep"] = ice_sep.collapsed('t',iris.analysis.MEAN)

# 3. SAT,d18O,precip 
variables=['d18O','temp','precip']
for ID in range(len(expID)):
	name = str(expName[ID])
	print ID, ":", name
	filename=''+data_root+expID[ID]+'.SAT_P_d18O.25kmx25km.nc'
	for VAR in variables:
		print "load vars:", str(VAR)
                if VAR=='temp':
                   con = iris.Constraint(cube_func=lambda cube: cube.var_name == VAR+'_mm_srf')
                   cube=iris.load_cube(filename, constraint=con)
                   cube=cube[-tlen:].collapsed(['time'],iris.analysis.MEAN)
                   cube.data=N.ma.masked_invalid(cube.data)
		   locals()[str(name)+"_"+str(VAR)] = cube
                if VAR=='precip':
                   con = iris.Constraint(cube_func=lambda cube: cube.var_name == VAR+'_mm_srf')
                   cube=iris.load_cube(filename, constraint=con)
                   cube=cube[-tlen:].collapsed(['time'],iris.analysis.MEAN)            
                   cube=cube*60*30*48*30.4375
                   cube.data=N.ma.masked_invalid(cube.data)
		   locals()[str(name)+"_"+str('P')] = cube
                else:
                   con = iris.Constraint(cube_func=lambda cube: cube.var_name == 'dO18')
                   cube=iris.load_cube(filename, constraint=con)
                   cube=cube[-tlen:].collapsed(['time'],iris.analysis.MEAN)    
                   cube.data=N.ma.masked_invalid(cube.data)
		   locals()[str(name)+"_"+str(VAR)] = cube

# 4. Orography
for ID in range(len(expID)):
    name = str(expName[ID])
    print ID, ":", name
    filename=''+data_root_max+'/netcdf/'+expID[ID]+'_d18O_sat_z.nc'
    con = iris.Constraint(cube_func=lambda cube: cube.var_name == 'orog')
    locals()[str(name)+"_"+str('orog')] = iris.load_cube(filename, constraint=con)

# 5. Mean Sea level pressure
for ID in range(len(expID)): 
    name =  str(expName[ID])
    print name
    print ID, ":", name
    root_to_data=''+data_root_max+'/'+expID[ID]+'/monthly'
    filename = ''+root_to_data+'/'+expID[ID]+'.p_mm_msl.monthly.nc'
    con = iris.Constraint(cube_func=lambda cube: cube.var_name == 'p_mm_msl')
    cube = iris.load_cube(filename, constraint=con)
    hp = cube[-tlen:].collapsed(['t','msl'],iris.analysis.MEAN)
    hp.data=N.ma.masked_invalid(hp.data)
    hp = hp /100
    locals()[str(name)+"_mslp"] = hp

# Defition of the reference simulation, i.e. for present days
reforog=PI_orog
refiso=PI_d18O
reftemp=PI_temp
refP=PI_P
refmslp=PI_mslp

# Arrays of simulations to plot
plotIDs=['BC128_1000lo','BC128_500lo','BC128','BC128_500hi','BC128_1000hi']
orogvars=[]
isovars=[]
satvars=[]
sivars=[]
Pvars=[]
mslpvars=[]

for ii in [expName[expName.index(x)] for x in plotIDs]:
    orogvars.append(locals()[str(ii)+"_orog"])
    isovars.append(locals()[str(ii)+"_d18O"])
    satvars.append(locals()[str(ii)+"_temp"])
    sivars.append(locals()[str(ii)+"_sep"])
    Pvars.append(locals()[str(ii)+"_P"])
    mslpvars.append(locals()[str(ii)+"_mslp"])


# ---------------------------------------------------------------
# COLORMAPS
# ---------------------------------------------------------------
import string
import matplotlib.colors as mcolors
colors1 = plt.cm.YlOrBr_r(np.linspace(0., 1, 128))
colors2 = plt.cm.BuPu(np.linspace(0, 1, 128))
colors = np.vstack((colors1, colors2))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
# colormap for absolute T
red = N.array([0, 0, 0, 10, 0, 0, 0, 60, 72, 127, 155, 221, 229, 239, 228, 205, 161, 116, 77]) / 256.
green = N.array([0, 0, 66, 144, 170, 191, 206, 209, 224, 255, 255, 242, 235, 190, 128, 72, 33, 29, 30]) / 256.
blue = N.array([128, 255, 255, 255, 255, 255, 209, 204, 208, 212, 250, 243, 99, 63, 39, 27, 22, 29, 27]) / 256.
colormap0 = N.array([red, green, blue]).T
my_colormap0 = matplotlib.colors.ListedColormap(colormap0, name='my_colormap0')
# colormap for T anomalies
red = N.array([0, 0, 0, 10, 0, 60, 72, 221, 255,255, 250, 255, 229, 239, 228, 205, 161, 116, 77]) / 256.
green = N.array([0, 0, 66, 144, 191, 209, 255, 242, 255,255, 250, 255, 235, 190, 128, 72, 33, 29, 30]) / 256.
blue = N.array([128, 255, 255, 255, 255, 204, 250, 243, 255, 255, 210, 0, 99, 63, 39, 27, 22, 29, 27]) / 256.
colormap = N.array([red, green, blue]).T
my_colormap = matplotlib.colors.ListedColormap(colormap, name='my_colormap')
fmt='%.0f'


# ---------------------------------------------------------------------------------------
# PLOT dZ, dP, dSAT, dd18O
# ---------------------------------------------------------------------------------------
filname='AIS_dZdPdSATd18Omslp_d1km_refPI_30yr_31Jan20' # Filename to save the plot

# Level for the sea-ice extent (10%)
icelevs = [0.15]

# Subplot indices for the different variables
idxz=[1,5,9,13,17]
idxP=[2,6,10,14,18]
idxtemp=[3,7,11,15,19]
idxiso=[4,8,12,16,20]

# Extraction of latitude and longitude
lat,lon=refmslp.coord('latitude').points,refmslp.coord('longitude').points


lon0 = 180
fig = plt.figure(figsize=(15,20))
Nsims=5
Nvars=4
sea_ice_red = ['7.6 %','1.1 %','','-6.0 %','-10.8 %'] #Sea-ice extent labels
rows = ['-1km', '-500m', 'no scaling', '+500m','+1km'] # Experiment labels (for each row)

#Orograpy subplots
for ID in range(len(isovars)):
	if ID == 2:
		dlevs=np.arange(0,4000.1,200)
		colormap = cm.viridis
		var1=(orogvars[ID].data)*lsm.data
		extent1=sivars[ID].data
		ax = fig.add_subplot(Nsims,Nvars,idxz[ID])
		var,lonsout = mpl_toolkits.basemap.shiftgrid(lon0, var1, lon, start=False, cyclic=360.0)
		var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
		extent,lonsout = mpl_toolkits.basemap.shiftgrid(lon0, extent1, lon, start=False, cyclic=360.0)
		extent,lonsout = mpl_toolkits.basemap.addcyclic(extent, lonsout)
		map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
		map.drawcoastlines()
		map.drawcountries()
		map.drawparallels(N.arange(-80.,81.,20.), color='grey')
		map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
		x,y = map(*np.meshgrid(lonsout,lat))
		contour_result5 = map.contourf(x, y, var, dlevs,cmap=colormap,extend='both')
		cax = map.contour(x,y,extent,icelevs,linestyles='-',colors='k',alpha=0.6)
		xpt,ypt = map(lon_ice,lat_ice)
		lonpt, latpt = map(xpt,ypt,inverse=True)
		collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
		ax.annotate(list(string.ascii_uppercase)[idxz[ID]-1], xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
		ax.annotate(sea_ice_red[ID],xy=(1,0.95),xycoords='axes fraction',color='grey',ha='right',va='top',fontsize='12',weight='bold')
                plt5_ax = plt.gca()
		colorbar_axes = fig.add_axes([0.08, 0.42, 0.005, 0.2])
		cbar = plt.colorbar(contour_result5, orientation= 'vertical', cax=colorbar_axes)
		cbar.set_label(r'z (km)',fontsize=18) # NOTE: NEED TO SWAP SIDE OF LABEL AND TICK LABELS ***
		cbar.set_ticks([0,2000,4000])
		cbar.set_ticklabels([0,2,4])
		cbar.ax.tick_params(length=0,labelsize=9)
		plt.gca().yaxis.set_ticks_position('left')
		plt.gca().yaxis.set_label_position('left')
	else:
		colormap = mymap
		dlevs=np.arange(-1000,1000.1,100)
		var1=(orogvars[ID].data-reforog.data)*lsm.data
		extent1=sivars[ID].data
		ax = fig.add_subplot(Nsims,Nvars,idxz[ID])
		var,lonsout = mpl_toolkits.basemap.shiftgrid(lon0, var1, lon, start=False, cyclic=360.0)
		var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
		extent,lonsout = mpl_toolkits.basemap.shiftgrid(lon0, extent1, lon, start=False, cyclic=360.0)
		extent,lonsout = mpl_toolkits.basemap.addcyclic(extent, lonsout)
		map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
		map.drawcoastlines()
		map.drawcountries()
		map.drawparallels(N.arange(-80.,81.,20.), color='grey')
		map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
		x,y = map(*np.meshgrid(lonsout,lat))
		contour_result1 = map.contourf(x, y, var, dlevs,cmap=colormap,extend='both')
		cax = map.contour(x,y,extent,icelevs,linestyles='-',colors='k',alpha=0.6)
		xpt,ypt = map(lon_ice,lat_ice)
		lonpt, latpt = map(xpt,ypt,inverse=True)
		collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
		ax.annotate(list(string.ascii_uppercase)[idxz[ID]-1], xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
                ax.annotate(sea_ice_red[ID],xy=(1,0.95),xycoords='axes fraction',color='grey',ha='right',va='top',fontsize='11',weight='bold')
		plt0_ax = plt.gca()

#dO18 subplots
n=0
dlevs=N.arange(-9,9.1,1.0)*1.4
dlevs2= N.arange(-10,10.1,2)
colormap = my_colormap
for ID in range(len(isovars)):
    var1=(isovars[ID].data-refiso.data)
    var1bis=(mslpvars[ID].data-refmslp.data)*lsm.data #mslp
    ax = fig.add_subplot(Nsims,Nvars,idxiso[ID])
    #var,lonsout = mpl_toolkits.basemap.shiftgrid(lon0, var1, lon_regridded, start=False, cyclic=360.0)
    #var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
    #for mslp
    var2,lonsout2 = mpl_toolkits.basemap.shiftgrid(lon0, var1bis, lon, start=False,cyclic=360.0)
    var2,lonsout2 = mpl_toolkits.basemap.addcyclic(var2, lonsout2)	
    map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
    map.drawparallels(N.arange(-80.,81.,20.), color='grey')
    map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
    map.drawcoastlines()
    map.drawcountries()
    x,y = map(lon_regridded,lat_regridded)
    x2,y2 = map(*np.meshgrid(lonsout2,lat))
    contour_result2 = map.contourf(x, y, var1[0,:,:], dlevs,cmap=colormap.reversed(),extend='both')
    cs = map.contour(x2,y2,var2,levels=dlevs2,colors='grey',linestyles='solid') #mslp contour
    xpt,ypt = map(lon_ice,lat_ice)
    lonpt, latpt = map(xpt,ypt,inverse=True)
    collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
    ax.annotate(list(string.ascii_uppercase)[idxiso[ID]-1], xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
    ax.set_ylabel(rows[n], size=16,rotation=0,labelpad=50)
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    plt1_ax = plt.gca()
    n=n+1

#SAT
for ID in range(len(isovars)):
    var1=(satvars[ID].data-reftemp.data)
    ax = fig.add_subplot(Nsims,Nvars,idxtemp[ID])
    #var,lonsout = mpl_toolkits.basemap.shiftgrid(lon0, var1, lon_regridded, start=False, cyclic=360.0)
    #var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)
    map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
    map.drawparallels(N.arange(-80.,81.,20.), color='grey')
    map.drawmeridians(N.arange(-180.,181.,20.), color='grey')
    map.drawcoastlines()
    map.drawcountries()
    x,y = map(lon_regridded,lat_regridded)
    contour_result3 = map.contourf(x, y, var1[0,:,:], dlevs,cmap=colormap,extend='both')
    xpt,ypt = map(lon_ice,lat_ice)
    lonpt, latpt = map(xpt,ypt,inverse=True)
    collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
    ax.annotate(list(string.ascii_uppercase)[idxtemp[ID]-1], xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
    plt2_ax = plt.gca()
        

#Precip & mean sea level pressure subplots
dlevs = N.arange(-10,10.1,0.1)
dlevs2= N.arange(-10,10.1,2)
for ID in range(len(isovars)):
    var1=(Pvars[ID].data-refP.data)*lsm_regridded.data
    #var1bis=(mslpvars[ID].data-refmslp.data)*lsm.data #mslp
    ax = fig.add_subplot(Nsims,Nvars,idxP[ID])
    #var,lonsout = mpl_toolkits.basemap.shiftgrid(lon0, var1, lon_regridded, start=False, cyclic=360.0)
    #var,lonsout = mpl_toolkits.basemap.addcyclic(var, lonsout)        
    #for mslp
    #var2,lonsout2 = mpl_toolkits.basemap.shiftgrid(lon0, var1bis, lon, start=False,cyclic=360.0)
    #var2,lonsout2 = mpl_toolkits.basemap.addcyclic(var2, lonsout2)	
    map = Basemap(projection='spaeqd',boundinglat=-55,lon_0=180,resolution='l')
    map.drawparallels(N.arange(-80.,81.,20.), color='grey')
    map.drawmeridians(N.arange(-180.,181.,20.), color='grey')    
    map.drawcoastlines()
    map.drawcountries()
    x,y = map(lon_regridded,lat_regridded)
    #x2,y2 = map(*np.meshgrid(lonsout2,lat))   
    contour_result4 = map.contourf(x, y, var1, dlevs,cmap=colormap.reversed(),extend='both')
    #cs = map.contour(x2,y2,var2,levels=dlevs2,colors='grey') #mslp contour   
    xpt,ypt = map(lon_ice,lat_ice)
    lonpt, latpt = map(xpt,ypt,inverse=True)
    collection = map.scatter(xpt,ypt, c='none',s=20, marker='o', edgecolors='k')
    ax.annotate(list(string.ascii_uppercase)[idxP[ID]-1], xy=(0, 1), xycoords='axes fraction',weight='bold',fontsize=14,xytext=(5, -5), textcoords='offset points',ha='left', va='top')
    plt3_ax = plt.gca()

left, bottom, width, height = plt0_ax.get_position().bounds
first_plot_left = plt0_ax.get_position().bounds[0]
width = left - first_plot_left + width
colorbar_axes = fig.add_axes([first_plot_left-0.02, bottom - 0.05, 0.1, 0.01])
cbar = plt.colorbar(contour_result1, colorbar_axes, orientation= 'horizontal')
cbar.set_label(r'$\Delta$z (m)',fontsize=14)
cbar.set_ticks([-1000,-400,0,400,1000])
cbar.set_ticklabels([-1000,-400,0,400,1000])
cbar.ax.tick_params(length=0,labelsize=12)

left, bottom, width, height = plt1_ax.get_position().bounds
first_plot_left = plt1_ax.get_position().bounds[0]
width = left - first_plot_left + width
colorbar_axes = fig.add_axes([first_plot_left+0.02, bottom - 0.05, 0.1, 0.01])
cbar = plt.colorbar(contour_result2, colorbar_axes, orientation= 'horizontal')
cbar.set_label(r'$\Delta\delta^{18}$O'+u' (\u2030)',fontsize=14)
cbar.set_ticks([-10,-5,0,5,10])
cbar.set_ticklabels([-10,-5,0,5,10])
cbar.ax.tick_params(length=0,labelsize=12)

left, bottom, width, height = plt2_ax.get_position().bounds
first_plot_left = plt2_ax.get_position().bounds[0]
width = left - first_plot_left + width
colorbar_axes = fig.add_axes([first_plot_left+0.005, bottom - 0.05, 0.1, 0.01])
cbar = plt.colorbar(contour_result3, colorbar_axes, orientation= 'horizontal')
cbar.set_label(r'$\Delta$SAT ($^{\circ}$C)',fontsize=14)
cbar.set_ticks([-10,-5,0,5,10])
cbar.set_ticklabels([-10,-5,0,5,10])
cbar.ax.tick_params(length=0,labelsize=12)

left, bottom, width, height = plt3_ax.get_position().bounds
first_plot_left = plt3_ax.get_position().bounds[0]
width = left - first_plot_left + width
colorbar_axes = fig.add_axes([first_plot_left-0.005, bottom -0.05, 0.1, 0.01])
cbar = plt.colorbar(contour_result4, colorbar_axes, orientation= 'horizontal')
cbar.set_label(r'$\Delta$P (mm m'+r'$^{-1}$)',fontsize=16)
cbar.set_ticks([-10,-5,0,5,10])
cbar.set_ticklabels([-10,-5,0,5,10])
cbar.ax.tick_params(length=0,labelsize=12)

plt.tight_layout()
fig.subplots_adjust(bottom=0.08)
fig.subplots_adjust(left=0.05)
plt.savefig(''+fig_root+'/'+filname+'4.pdf', dpi=300)
plt.savefig(''+fig_root+'/'+filname+'4.png', dpi=300)

#End
