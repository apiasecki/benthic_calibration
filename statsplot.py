# This takes an excell plot that has the bootsrap statistics output from
# the bjarte script and just replots them in python for prettiness
# the excell file has two sheets for every species
# sheet 1 = species specific x, y and fit
# sheet 2 = all species x, y, fit, 0.05%, 0.95%

import pandas
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *
from pyteomics import mass # can do cool isotope math stuff, not used here
# plt.style.use('ggplot')
#plt.style.use('presentation')     #have this when i want slides, basically makes it readable from the back row

# import file and make main variables
xls_file = pandas.ExcelFile('stats_02_13_18.xlsx')
all1  = xls_file.parse('All Data1')
all2  = xls_file.parse('All Data2')
cibs1 = xls_file.parse('Cibicidoides1')
cibs2 = xls_file.parse('Cibicidoides2')
hoe1  = xls_file.parse('Hoeglundina1')
hoe2  = xls_file.parse('Hoeglundina2')
len1  = xls_file.parse('Lenticulina1')
len2  = xls_file.parse('Lenticulina2')
pyr1  = xls_file.parse('Pyrgo1')
pyr2  = xls_file.parse('Pyrgo2')
uvi1  = xls_file.parse('Uvigerina1')
uvi2  = xls_file.parse('Uvigerina2')


# plot it up
fig, ((ax1, ax5), (ax2, ax4), (ax6, ax3))= plt.subplots(3,2, sharex='col', sharey='row', figsize = (8,10))
plt.rc('font', family = 'Helvetica')
# ax1 = fig.add_subplot(1, 3, 2)
ax1.plot(all1['xtax'], all1['ytax'], '+', color = 'k', label = 'All Data')
ax1.plot(all1['xtax'], all1['txline'], '-', color = 'k')
ax1.plot(all2['x'], all2['low'], '-.', color = 'k')
ax1.plot(all2['x'], all2['high'], '-.', color = 'k')
# ax1.legend(loc = 'lower right')
ax1.set_title('All Data')
# ax1.set_xlabel('$10^6/$T$^2 ($K$)$')
ax1.set_ylabel('$\Delta_{47}$ (\u2030)')

# ax2 = fig.add_subplot(2, 3, 2)
ax2.plot(cibs1['xtax'], cibs1['ytax'], 'o', color = 'r', label = 'C. pachyderma')
ax2.plot(cibs1['xtax'], cibs1['txline'], '-', color = 'r', label = 'C. pachyderma fit')
ax2.plot(cibs2['x'], cibs2['y'], '+', color = '0.75', zorder = 1)
ax2.plot(cibs2['x'], cibs2['median'], '-', color = 'k', label = 'All Data fit')
ax2.plot(cibs2['x'], cibs2['low'], '-.', color = 'k')
ax2.plot(cibs2['x'], cibs2['high'], '-.', color = 'k')
# ax2.legend(loc = 'lower right')
ax2.set_title('C. pachyderma')
# ax2.set_xlabel('$10^6/$T$^2 ($K$)$')
ax2.set_ylabel('$\Delta_{47}$ (\u2030)')

# ax3 = fig.add_subplot(3, 3, 2)
ax3.plot(hoe1['xtax'], hoe1['ytax'], 'p', color = 'darkviolet', label = 'H. elegans')
ax3.plot(hoe1['xtax'], hoe1['txline'], '-', color = 'darkviolet', label = 'H. elegans fit')
ax3.plot(hoe2['x'], hoe2['y'], '+', color = '0.75', zorder = 1)
ax3.plot(hoe2['x'], hoe2['median'], '-', color = 'k', label = 'All Data fit')
ax3.plot(hoe2['x'], hoe2['low'], '-.', color = 'k')
ax3.plot(hoe2['x'], hoe2['high'], '-.', color = 'k')
# ax3.legend(loc = 'lower right')
ax3.set_title('H. elegans')
ax3.set_xlabel('$10^6/$T$^2 ($K$)$')
# ax3.set_ylabel('$\Delta_{47}$ (\u2030)')

# ax4 = fig.add_subplot(4, 3, 2)
ax4.plot(len1['xtax'], len1['ytax'], '^', color = 'g', label = 'Lenticulina spp')
ax4.plot(len1['xtax'], len1['txline'], '-', color = 'g', label = 'Lenticulina spp fit')
ax4.plot(len2['x'], len2['y'], '+', color = '0.75', zorder = 1)
ax4.plot(len2['x'], len2['median'], '-', color = 'k', label = 'All Data fit')
ax4.plot(len2['x'], len2['low'], '-.', color = 'k')
ax4.plot(len2['x'], len2['high'], '-.', color = 'k')
# ax4.legend(loc = 'lower right')
ax4.set_title('Lenticulina spp')
# ax4.set_xlabel('$10^6/$T$^2 ($K$)$')
# ax4.set_ylabel('$\Delta_{47}$ (\u2030)')

# ax5 = fig.add_subplot(5, 3, 2)
ax5.plot(pyr1['xtax'], pyr1['ytax'], 'h', color = 'c', label = 'Pyrgo spp')
ax5.plot(pyr1['xtax'], pyr1['txline'], '-', color = 'c', label = 'Pyrgo spp fit')
ax5.plot(pyr2['x'], pyr2['y'], '+', color = '0.75', zorder = 1)
ax5.plot(pyr2['x'], pyr2['median'], '-', color = 'k', label = 'All Data fit')
ax5.plot(pyr2['x'], pyr2['low'], '-.', color = 'k')
ax5.plot(pyr2['x'], pyr2['high'], '-.', color = 'k')
# ax5.legend(loc = 'lower right')
ax5.set_title('Pyrgo spp')
# ax5.set_xlabel('$10^6/$T$^2 ($K$)$')
# ax5.set_ylabel('$\Delta_{47}$ (\u2030)')

# ax6 = fig.add_subplot(6, 3, 2)
ax6.plot(uvi1['xtax'], uvi1['ytax'], 's', color = 'b', label = 'U. mediteranea')
ax6.plot(uvi1['xtax'], uvi1['txline'], '-', color = 'b', label = 'U. mediteranea')
ax6.plot(uvi2['x'], uvi2['y'], '+', color = '0.75', zorder = 1)
ax6.plot(uvi2['x'], uvi2['median'], '-', color = 'k', label = 'All Data fit')
ax6.plot(uvi2['x'], uvi2['low'], '-.', color = 'k')
ax6.plot(uvi2['x'], uvi2['high'], '-.', color = 'k')
# ax6.legend(loc = 'lower right')
ax6.set_title('U. mediteranea')
ax6.set_xlabel('$10^6/$T$^2 ($K$)$')
ax6.set_ylabel('$\Delta_{47}$ (\u2030)')

plt.savefig('02_15_18_statsall.pdf')
