import matplotlib.pyplot as plt
import numpy as np

convert_GeV_pbarn = (2.57*10**-9)


#readout data from files
s12_ana = np.fromfile('sigma_ana_s12.bytes', dtype=np.float64)
sigma_ana = np.fromfile('sigma_ana.bytes', dtype=np.float64)

s12_num = np.fromfile('sigma_num_s12.bytes', dtype=np.float64)
sigma_num = np.fromfile('sigma_num.bytes', dtype=np.float64)

s12_hT = np.fromfile('sigma_num_hT_s12.bytes', dtype=np.float64)
sigma_hT = np.fromfile('sigma_num_hT.bytes', dtype=np.float64)

s12_ana_hT = np.fromfile('sigma_ana_hT_s12.bytes', dtype=np.float64)
sigma_ana_hT = np.fromfile('sigma_ana_hT.bytes', dtype=np.float64)


#make plots
plt.plot(np.sqrt(s12_ana), sigma_ana/convert_GeV_pbarn, label=r'$\sigma_{\mathrm{ana}}$', color='darkgreen')
plt.plot(np.sqrt(s12_num), sigma_num/convert_GeV_pbarn, label=r'$\sigma_{\mathrm{num}}$', linestyle='dashed', color='lightgreen')
plt.plot(np.sqrt(s12_hT), sigma_hT/convert_GeV_pbarn, label=r'$\sigma_\mathrm{hT}$', color='blue')
plt.plot(np.sqrt(s12_ana_hT), sigma_ana_hT/convert_GeV_pbarn, label=r'$\sigma_\mathrm{expBReg\_hT}$', color='lightblue', linestyle='dashed')

#formatting plot
plt.xlabel(r'$\sqrt{s_{12}}$ in GeV')
plt.ylabel(r'$\mathrm{d}\sigma/\mathrm{d}s_{12}$ in pbarn')
plt.yscale('log')
plt.yticks([1, 2, 5, 10, 20, 50, 100], [1, 2, 5, 10, 20, 50, 100])
plt.legend()
plt.savefig('sigma.pdf')
