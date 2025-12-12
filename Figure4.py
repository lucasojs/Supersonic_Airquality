import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import math

ds_area = xr.open_dataset("/home/jsoh/regular_lat_lon_91x144.nc")
AREA_91_144= np.array(ds_area['area'][:])
ds_pop = xr.open_dataset("/home/jsoh/Population_2015_2x2.5_PC_DC_global.nc")
Pop_91_144 = np.array(ds_pop['pop'][0,:,:])

ds = xr.open_dataset("/home/jsoh/Air_quality_paper/surface_fields_by_cases.nc")
ds2 =  xr.open_dataset("/home/jsoh/Air_quality_paper/surface_fields_base.nc")


x = np.array(ds['lon'][:])
y = np.array(ds['lat'][:])
X, Y = np.meshgrid(y, x)

PM_25_total_cases= ds['sfc_PM25'].values
PM_25_N_total_cases=ds['sfc_PM25_nitrate'].values
PM_25_S_total_cases=ds['sfc_PM25_sulfate'].values
PM_25_A_total_cases=ds['sfc_PM25_ammo'].values

PM_25_N_base = ds2['sfc_PM25_nitrate'].values
PM_25_S_base = ds2['sfc_PM25_sulfate'].values
PM_25_A_base = ds2['sfc_PM25_ammo'].values
PM_25_base = ds2['sfc_PM25'].values

PWM_PM = np.zeros((35))
PWM_PM_N = np.zeros((35))
PWM_PM_S = np.zeros((35))
PWM_PM_A = np.zeros((35))

for i in range(35):
    PWM_PM[i] = np.sum((PM_25_total_cases[i, :, :]-PM_25_base) * Pop_91_144) / np.sum(Pop_91_144)
    PWM_PM_N[i] = np.sum((PM_25_N_total_cases[i, :, :]-PM_25_N_base) * Pop_91_144) / np.sum(Pop_91_144)
    PWM_PM_S[i] = np.sum((PM_25_S_total_cases[i, :, :]-PM_25_S_base) * Pop_91_144) / np.sum(Pop_91_144)
    PWM_PM_A[i] = np.sum((PM_25_A_total_cases[i, :, :]-PM_25_A_base) * Pop_91_144) / np.sum(Pop_91_144)

AVG_PM = np.zeros((35))
AVG_PM_N = np.zeros((35))
AVG_PM_S = np.zeros((35))
AVG_PM_A = np.zeros((35))

for i in range(35):
    AVG_PM[i] = np.sum((PM_25_total_cases[i, :, :]-PM_25_base) * AREA_91_144) / np.sum(AREA_91_144)
    AVG_PM_N[i] = np.sum((PM_25_N_total_cases[i, :, :]-PM_25_N_base) * AREA_91_144) / np.sum(AREA_91_144)
    AVG_PM_S[i] = np.sum((PM_25_S_total_cases[i, :, :]-PM_25_S_base) * AREA_91_144) / np.sum(AREA_91_144)
    AVG_PM_A[i] = np.sum((PM_25_A_total_cases[i, :, :]-PM_25_A_base) * AREA_91_144) / np.sum(AREA_91_144)
    
portion_N = PWM_PM_N / 3
portion_S = PWM_PM_S / 3
portion_A = PWM_PM_A / 3
portion_T = PWM_PM / 3

max_abs_diff = 0.85
selected_indices = [1, 6, 11, 16, 21, 26, 31]
altitude_labels = ['8–10', '10–12', '12–14', '14–16', '16–18', '18–20', '20–22']
index = np.arange(len(selected_indices))

portion_N_selected = portion_N[selected_indices]
portion_S_selected = portion_S[selected_indices]
portion_A_selected = portion_A[selected_indices]
portion_T_selected = portion_T[selected_indices]

fig, ax = plt.subplots(figsize=(2, 1.55))

ax.plot(portion_N_selected, index, '-.', label='Nitrate', color='darkred', linewidth=1)
ax.plot(portion_S_selected, index, '--', label='Sulfate', color='darkblue', linewidth=1)
ax.plot(portion_A_selected, index, ':', label='Ammonium', color='orange', linewidth=1)
ax.plot(portion_T_selected, index, '-', label='PM₂.₅', color='black', linewidth=1)
ax.set_yticks(index)
ax.set_yticklabels(altitude_labels, fontsize=6)
ax.tick_params(axis='x', labelsize=6)
ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
ax.legend(fontsize=6, loc='lower right', frameon=False)

plt.tight_layout()
plt.savefig('./Figure4_Submission/Figure4_a_only.png', dpi=600, bbox_inches='tight')

fig = plt.figure(figsize=(3, 1.55))  # 약 23cm x 16.5cm
gs = gridspec.GridSpec(nrows=2, ncols=3, width_ratios=[1, 1, 0.05], height_ratios=[1, 1], hspace=0.02, wspace=0.01)

index1 = 1
index2 = 31

cases = [
    (PM_25_N_total_cases[index1], PM_25_N_base),
    (PM_25_N_total_cases[index2], PM_25_N_base),
    (PM_25_S_total_cases[index1], PM_25_S_base),
    (PM_25_S_total_cases[index2], PM_25_S_base)
]

for i, (case_data, base_data) in enumerate(cases):
    row = i // 2
    col = i % 2
    ax = fig.add_subplot(gs[row, col], projection = ccrs.Robinson())

    diff = (case_data - base_data) / 3
    im = ax.pcolormesh(
        Y, X, np.transpose(diff),
        transform=ccrs.PlateCarree(),
        vmin=-max_abs_diff, vmax=max_abs_diff,
        cmap='RdBu_r'
    )
    ax.add_feature(cfeature.COASTLINE.with_scale('110m'))
    
cbar_ax = fig.add_axes([0.89, 0.11, 0.02, 0.78])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
cbar.ax.tick_params(labelsize=8)

plt.savefig('./Figure4_Submission/Figure4_BCDE_only.png', dpi=600, bbox_inches='tight')