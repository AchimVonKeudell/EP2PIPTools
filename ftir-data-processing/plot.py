import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# fn = r"D:\Promotie\data\output\23_IRRAS-NOx\20250319_IRRAS-NOx-vary-ratio\FTIR\5_N2+O2_plasma\NO\19sccmN2+1sccmO2--0000.dat"
fn = r"D:\Promotie\data\output\20251216_IRRAS-N2+H2+CH4\FTIR\2_N2+H2+CH4_plasma\HCN\20sccmN2+1sccmH2+1sccmCH4.0024.dat"
data = np.loadtxt(fn)

fig = plt.figure(figsize=(5, 3))
# gs = gridspec.GridSpec(3, 1)
ax = fig.add_subplot()
# ax.set(xlabel='wavenumber / cm$^{-1}$', xlim=(2800, 3200),
#        ylabel='transmittance / -')
# ax1 = fig.add_subplot(gs[2], sharex=ax)
# ax1.set(, ylabel='residual / -')

ax.plot(data[:, 0], data[:, 1], label='measured')
ax.plot(data[:, 0], data[:, 2], label='fit')
ax.plot(data[:, 0], data[:, 1] - data[:, 2] + 1, c='black', label='residual')

ax.legend(loc=3)
# ax.plot(data[:, 0], -1e3 * np.log(data[:, 1]))
fig.tight_layout()
# fig.savefig(r"D:\Promotie\data\output\20251120_IRRAS-N2+CH4-diluted\FTIR\3_Ar+CH4_plasma\exmaple_CH4.png", dpi=200)
plt.show()