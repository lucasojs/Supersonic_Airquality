import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

ds = xr.open_dataset("/home/jsoh/Air_quality_paper/surface_fields_by_cases.nc")
z_total_cases = ds['sfc_O3'].values

ds2 =  xr.open_dataset("/home/jsoh/Air_quality_paper/surface_fields_base.nc")
O3_base = ds2['sfc_O3'].values

ds3 = xr.open_dataset("/home/jsoh/Air_quality_paper/Fields_3cases.nc")

diff_JO3_case17 = ds3['diff_JO3'].values[1,:]
diff_JO3_case32 = ds3['diff_JO3'].values[2,:]

diff_O3_case17 = (z_total_cases[16,:]-O3_base)/3
diff_O3_case32 = (z_total_cases[31,:]-O3_base)/3

O3_base = O3_base

diff_O3_flat_17 = diff_O3_case17.flatten()
JO3_O3_flat_17 = (O3_base * diff_JO3_case17).flatten()

diff_O3_flat_32 = diff_O3_case32.flatten()
JO3_O3_flat_32 = (O3_base * diff_JO3_case32).flatten()

fig, ax = plt.subplots(1, 2, figsize=(3.54, 1.8))
plt.subplots_adjust(wspace=0.3)

hb1 = ax[0].hexbin(JO3_O3_flat_17 * 1e3, diff_O3_flat_17, gridsize=35, cmap='plasma', mincnt=10, linewidths=0)
ax[0].axhline(0, color='black', linestyle='dashed', linewidth=0.8)
ax[0].axvline(0, color='black', linestyle='dashed', linewidth=0.8)
ax[0].tick_params(labelsize=8)

hb2 = ax[1].hexbin(JO3_O3_flat_32 * 1e3, diff_O3_flat_32, gridsize=35, cmap='plasma', mincnt=10, linewidths=0)
ax[1].axhline(0, color='black', linestyle='dashed', linewidth=0.8)
ax[1].axvline(0, color='black', linestyle='dashed', linewidth=0.8)
ax[1].tick_params(labelsize=8)

cbar = fig.colorbar(hb2, ax=ax, orientation='vertical', fraction=0.025, pad=0.02)
cbar.ax.tick_params(labelsize=8)

plt.savefig('./Figure3_Submission/Figure3.png', dpi=600, bbox_inches='tight')