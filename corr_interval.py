import pandas
import numpy as np
import matplotlib.pyplot as plt

xls_file = pandas.ExcelFile('corr_interval.xlsx')
data = xls_file.parse('Sheet1')


fig1 = plt.figure()
plt.rc('font', family='Helvetica')
ax = fig1.add_subplot(1,1,1)

ax.plot(data.num, data.IsoA, 'o-.', color = 'red', label = 'ETH1')
ax.plot(data.num, data.IsoB, 'o-.', color = 'blue', label = 'ETH2 (tracer)')
ax.plot(data.num, data.Chalk, 'o-.', color = 'pink', label = 'ETH3')
ax.plot(data.num, data.Riedel, 'o-.', color = 'green', label = 'ETH4')
plt.axvline(x = 40, color = 'gray', linewidth= 0.5, linestyle = ':')

ax.set_xlabel('Number of standards used for correction')
ax.set_ylabel('Standard Deviation $\Delta_{47}$ (\u2030)')
ax.legend(loc='best')

plt.savefig('corr_interval.pdf')
