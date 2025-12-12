import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import math
ds_area = xr.open_dataset("/home/jsoh/regular_lat_lon_91x144.nc")
AREA_91_144= np.array(ds_area['area'][:])
ds_pop = xr.open_dataset("/home/jsoh/Population_2015_2x2.5_PC_DC_global.nc")
Pop_91_144 = np.array(ds_pop['pop'][0,:,:])

ds = xr.open_dataset("/home/jsoh/Air_quality_paper/surface_fields_by_cases.nc")
z_total_cases = ds['sfc_O3'].values
PM_25_total_cases = ds['sfc_PM25'].values
ds2 =  xr.open_dataset("/home/jsoh/Air_quality_paper/Fields_base.nc")
z_base = ds2['sfc_O3'].values
PM_25_base = ds2['sfc_PM25'].values

index1 = 1
index3 = 31

x = np.array(ds['lon'][:])
y = np.array(ds['lat'][:])
X, Y = np.meshgrid(y, x)

diff1 = (z_total_cases[index1, :, :] - z_base) / 3
diff3 = (z_total_cases[index3, :, :] - z_base) / 3 

plt.rcParams.update({'font.size': 11})

fig = plt.figure(figsize=(5, 2)) 
gs = fig.add_gridspec(2, 2, height_ratios=[15, 1], hspace=0.01, wspace=0.05)

# (a)
ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.Robinson())
im1 = ax1.pcolormesh(Y, X, np.transpose(diff1), transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=-3.5, vmax=3.5)
ax1.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax1.set_global()

# (b)
ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.Robinson())
im2 = ax2.pcolormesh(Y, X, np.transpose(diff3), transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=-3.5, vmax=3.5)
ax2.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax2.set_global()

cax = fig.add_subplot(gs[1, :])
cbar = fig.colorbar(im2, cax=cax, orientation='horizontal',fraction=0.046, pad=0.07)
cbar.ax.tick_params(labelsize=8)

plt.savefig('./Figure1_Submission/Figure1_ab_only.png', dpi=600, bbox_inches='tight')

PWM_z = np.zeros((35))


for i in range(35):
    PWM_z[i] = np.sum((z_total_cases[i, :, :]-z_base) * Pop_91_144) / np.sum(Pop_91_144)
    

AVG_z = np.zeros((35))


for i in range(35):
    AVG_z[i] = np.sum((z_total_cases[i, :, :]-z_base)  * AREA_91_144) / np.sum(AREA_91_144)
    
    
sets_avg ={
    "NH": (AVG_z[[0, 5, 10, 15, 20, 25, 30]]*(1-np.sqrt(3)/2) + AVG_z[[1, 6, 11, 16, 21, 26, 31]]*(np.sqrt(3)/2-0.5) + AVG_z[[2, 7, 12, 17, 22, 27, 32]]*0.5)/3,
    "SH": (AVG_z[[3, 8, 13, 18, 23, 28, 33]]*0.5 + AVG_z[[4, 9, 14, 19, 24, 29, 34]]*0.5) / 3
}

sets_pwm = {
    "NH": (PWM_z[[0, 5, 10, 15, 20, 25, 30]]*(1-np.sqrt(3)/2) + PWM_z[[1, 6, 11, 16, 21, 26, 31]]*(np.sqrt(3)/2-0.5) + PWM_z[[2, 7, 12, 17, 22, 27, 32]]*0.5) / 3,
    "SH": (PWM_z[[3, 8, 13, 18, 23, 28, 33]]*0.5 + PWM_z[[4, 9, 14, 19, 24, 29, 34]]*0.5) / 3 
}

color_nh = "royalblue"  
color_sh = "darkorange" 

fig, ax = plt.subplots(figsize=(3, 2))

ax.plot(sets_avg["NH"], range(1, 8), label="AWM NH", linewidth=1, color=color_nh, linestyle='--')
ax.plot(sets_avg["SH"], range(1, 8), label="AWM SH", linewidth=1, color=color_sh, linestyle='--')
ax.plot(sets_pwm["NH"], range(1, 8), label="PWM NH", linewidth=1, color=color_nh)
ax.plot(sets_pwm["SH"], range(1, 8), label="PWM SH", linewidth=1, color=color_sh)

ax.set_yticks(range(1, 8))
ax.set_yticklabels(['8-10', '10-12', '12-14', '14-16', '16-18', '18-20', '20-22'], fontsize=8)
ax.set_ylim(0.5, 7.5)

ax.tick_params(axis='x', labelsize=8)
ax.set_xticks([-1.5,-1,-0.5,0,0.5])
ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)
ax.legend(fontsize=8, loc='lower left', frameon=False)

plt.tight_layout()
plt.savefig('./Figure1_Submission/Figure1_c_only.png', dpi=600, bbox_inches='tight')
