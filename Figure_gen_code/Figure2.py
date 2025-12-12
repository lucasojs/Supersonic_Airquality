import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ds = xr.open_dataset("/home/jsoh/Air_quality_paper/Fields_base.nc")
ds2 = xr.open_dataset("/home/jsoh/Air_quality_paper/Fields_3cases.nc")

x = np.array(ds['lon'][:])
y = np.array(ds['lat'][:])
height = np.array(ds['height'][:])

X, Y = np.meshgrid(y, x)
H, L = np.meshgrid(height, y)

N_lat_hbase = ds['zonal_O3'].values
N_lat_h_case2 = ds2['zonal_O3'].values[0,:]
N_lat_h_case17 = ds2['zonal_O3'].values[1,:]
N_lat_h_case32 = ds2['zonal_O3'].values[2,:]

diff2 = np.transpose(N_lat_h_case2 - N_lat_hbase) / np.transpose(N_lat_hbase) * 100 / 3
diff17 = np.transpose(N_lat_h_case17 - N_lat_hbase) / np.transpose(N_lat_hbase) * 100 / 3
diff32 = np.transpose(N_lat_h_case32 - N_lat_hbase) / np.transpose(N_lat_hbase) * 100 / 3

cmaps = ['RdBu_r', 'RdBu_r', 'PuOr_r']
diffs = [diff2, diff17, diff32]
L_cuts = [55, 55, 55]

fig = plt.figure(figsize=(7, 1.5))
gs = gridspec.GridSpec(nrows=2, ncols=3, height_ratios=[18, 1], width_ratios=[1, 1, 1.1], hspace=0.05, wspace=0.2)

ax1 = fig.add_subplot(gs[0, 0])
im1 = ax1.pcolormesh(L[:, :L_cuts[0]], H[:, :L_cuts[0]], diffs[0][:, :L_cuts[0]], vmin=-6, vmax=6, cmap=cmaps[0], shading='auto')
ax1.set_xticks([-90, -45,  0, 45, 90])
ax1.set_xticklabels([r'90°S', r'45°S',  r'Eq', r'45°N', r'90°N'], fontsize=8)
ax1.tick_params(axis='x', labelsize=8)
ax1.tick_params(axis='y', labelsize=8)
ax1.set_yticks([0, 10, 20, 30, 40])
ax1.set_yticklabels(['0', '10', '20', '30', '40'], fontsize=8)

ax2 = fig.add_subplot(gs[0, 1])
im2 = ax2.pcolormesh(L[:, :L_cuts[1]], H[:, :L_cuts[1]], diffs[1][:, :L_cuts[1]], vmin=-6, vmax=6, cmap=cmaps[1], shading='auto')
ax2.set_xticks([-90, -45,  0, 45, 90])
ax2.set_xticklabels([r'90°S', r'45°S',  r'Eq', r'45°N', r'90°N'], fontsize=8)
ax2.tick_params(axis='x', labelsize=8)
ax2.tick_params(axis='y', labelsize=8)
ax2.set_yticks([0, 10, 20, 30, 40])
ax2.set_yticklabels([])

cbar_ax12 = fig.add_axes([0.122, -0.05, 0.492, 0.05])
cbar_ab = fig.colorbar(im2, cax=cbar_ax12, orientation='horizontal')
cbar_ab.ax.tick_params(labelsize=8)

ax3 = fig.add_subplot(gs[0, 2])
im3 = ax3.pcolormesh(L[:, :L_cuts[2]], H[:, :L_cuts[2]], diffs[2][:, :L_cuts[2]], vmin=-40, vmax=40, cmap=cmaps[2], shading='auto')

ax3.set_xticks([-90, -45,  0, 45, 90])
ax3.set_xticklabels([r'90°S', r'45°S',  r'Eq', r'45°N', r'90°N'], fontsize=8)
ax3.tick_params(axis='x', labelsize=8)
ax3.tick_params(axis='y', labelsize=8)
ax3.set_yticklabels([])
ax1.set_yticklabels(['0', '10', '20', '30', '40'], fontsize=8)
ax3.set_yticklabels([])
ax3.set_yticks([0, 10, 20, 30, 40])

cbar_ax3 = fig.add_axes([0.66, -0.05, 0.24, 0.05])  
cbar_c = fig.colorbar(im3, cax=cbar_ax3, orientation='horizontal')
cbar_c.ax.tick_params(labelsize=8)

plt.savefig('./Figure2_Submission/Figure2_abc_only.png', dpi=600, bbox_inches='tight')

column_base = ds['col_O3'].values
column_case2 = ds2['col_O3'].values[0,:]
column_case17 = ds2['col_O3'].values[1,:]
column_case32 = ds2['col_O3'].values[2,:]

diff1 = (column_case2 - column_base) / 3 / column_base * 100
diff2 = (column_case17 - column_base) / 3 / column_base * 100
diff3 = (column_case32 - column_base) / 3 / column_base * 100
diffs = [np.transpose(diff1), np.transpose(diff2), np.transpose(diff3)]

vmins = [-2.2, -2.2, -20]
vmaxs = [2.2, 2.2, 20]

fig = plt.figure(figsize=(7, 1.5))
gs = gridspec.GridSpec(nrows=2, ncols=3, height_ratios=[18, 1], width_ratios=[1, 1, 1.1], hspace=0.05, wspace=0.2)

ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.Robinson())
im1 = ax1.pcolormesh(Y, X, diffs[0], transform=ccrs.PlateCarree(), cmap=cmaps[0], vmin=vmins[0], vmax=vmaxs[0])
ax1.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax1.set_global()


ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.Robinson())
im2 = ax2.pcolormesh(Y, X, diffs[1], transform=ccrs.PlateCarree(), cmap=cmaps[1], vmin=vmins[1], vmax=vmaxs[1])
ax2.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax2.set_global()

cbar_ax12 = fig.add_axes([0.122, 0.13, 0.492, 0.05])  # left, bottom, width, height
cbar_ab = fig.colorbar(im2, cax=cbar_ax12, orientation='horizontal')
cbar_ab.ax.tick_params(labelsize=8)

ax3 = fig.add_subplot(gs[0, 2], projection=ccrs.Robinson())
im3 = ax3.pcolormesh(Y, X, diffs[2], transform=ccrs.PlateCarree(), cmap=cmaps[2], vmin=vmins[2], vmax=vmaxs[2])
ax3.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax3.set_global()

cbar_ax3 = fig.add_axes([0.66, 0.13, 0.24, 0.05])
cbar_c = fig.colorbar(im3, cax=cbar_ax3, orientation='horizontal')
cbar_c.ax.tick_params(labelsize=8)

plt.savefig('./Figure2_Submission/Figure2_def_only.png', dpi=600, bbox_inches='tight')


diff_JO3_case2 = ds2['diff_JO3'].values[0,:]
diff_JO3_case17 = ds2['diff_JO3'].values[1,:]
diff_JO3_case32 = ds2['diff_JO3'].values[2,:]

diff1 = diff_JO3_case2*10**5
diff2 = diff_JO3_case17*10**5
diff3 = diff_JO3_case32*10**5
diffs = [np.transpose(diff1), np.transpose(diff2), np.transpose(diff3)]

cmaps = ['PRGn', 'PRGn', 'plasma']
vmins = [-0.18, -0.18, 0]
vmaxs = [0.18, 0.18, 3]


fig = plt.figure(figsize=(7, 1.5))
gs = gridspec.GridSpec(nrows=2, ncols=3, height_ratios=[18, 1], width_ratios=[1, 1, 1.1], hspace=0.05, wspace=0.2)

ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.Robinson())
im1 = ax1.pcolormesh(Y, X, diffs[0], transform=ccrs.PlateCarree(), cmap=cmaps[0], vmin=vmins[0], vmax=vmaxs[0])
ax1.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax1.set_global()


ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.Robinson())
im2 = ax2.pcolormesh(Y, X, diffs[1], transform=ccrs.PlateCarree(), cmap=cmaps[1], vmin=vmins[1], vmax=vmaxs[1])
ax2.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax2.set_global()

cbar_ax12 = fig.add_axes([0.122, 0.13, 0.492, 0.05])  # left, bottom, width, height
cbar_ab = fig.colorbar(im2, cax=cbar_ax12, orientation='horizontal')
cbar_ab.ax.tick_params(labelsize=8)

ax3 = fig.add_subplot(gs[0, 2], projection=ccrs.Robinson())
im3 = ax3.pcolormesh(Y, X, diffs[2], transform=ccrs.PlateCarree(), cmap=cmaps[2], vmin=vmins[2], vmax=vmaxs[2])
ax3.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax3.set_global()

cbar_ax3 = fig.add_axes([0.66, 0.13, 0.24, 0.05])
cbar_c = fig.colorbar(im3, cax=cbar_ax3, orientation='horizontal')
cbar_c.ax.tick_params(labelsize=8)

plt.savefig('./Figure2_Submission/Figure2_ghi_only.png', dpi=600, bbox_inches='tight')