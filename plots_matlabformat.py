# THIS TAKES AN ALREADY FORMATED DATA TABLE from matlab AND DOES THE REGRESSIONS
# IT ALSO MAKES THE PLOTS
# LAST EDITED 11-29-17

import pandas
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *
from pyteomics import mass # can do cool isotope math stuff, not used here
# plt.style.use('ggplot')
#plt.style.use('presentation')     #have this when i want slides, basically makes it readable from the back row

# import file and make main variables
xls_file = pandas.ExcelFile('out_02_13_18.xlsx')
avg = xls_file.parse('avg')
avg10 = xls_file.parse('avg10')
reps = xls_file.parse('all')
cibs = xls_file.parse('cibs')
lent = xls_file.parse('len')
pyrgo = xls_file.parse('pyrgo')
site = xls_file.parse('sites')
uvi = xls_file.parse('uvi')
elegans = xls_file.parse('helegans')
other = xls_file.parse('assorted')
# big = xls_file.parse('big')
# notbig = xls_file.parse('not')

otherpeople = pandas. ExcelFile('others.xlsx')
travertine = otherpeople.parse('travertine')
leeds = otherpeople.parse('leeds')
ETH = otherpeople.parse('ETH')
cavepearls = otherpeople.parse('cavepearls')
magali = otherpeople.parse('bonifacie')
#others.to_csv('others.csv')


# create 1/T^2 variable
# avg['XT'] = 1e6 /((273.15+avg['temp'])**2)
# avg10['XT'] = 1e6 /((273.15+avg10['temp'])**2)
# reps['XT'] = 1e6 /((273.15+reps['temp'])**2)
# cibs['XT'] = 1e6 /((273.15+cibs['temp'])**2)
# lent['XT'] = 1e6 /((273.15+lent['temp'])**2)
# pyrgo['XT'] = 1e6 /((273.15+pyrgo['temp'])**2)
# uvi['XT'] = 1e6 /((273.15+uvi['temp'])**2)
# elegans['XT'] = 1e6 /((273.15+elegans['temp'])**2)
# other['XT'] = 1e6 /((273.15+other['temp'])**2)
# big['XT'] = 1e6 /((273.15+big['temp'])**2)
# notbig['XT'] = 1e6 /((273.15+notbig['temp'])**2)
# site['XT'] = 1e6 /((273.15+site['temp'])**2)
travertine['XT'] = 1e6 /((273.15+travertine['temp'])**2)
leeds['XT'] = 1e6 /((273.15+leeds['temp'])**2)
ETH['XT'] = 1e6 /((273.15+ETH['temp'])**2)
cavepearls['XT'] = 1e6 /((273.15+cavepearls['temp'])**2)


# regression
def line(x, a, b):
    return a*x+b
popt_all, pcov_all = curve_fit(line, reps['XT'], reps['D47'])
#popt5, pcov5 = curve_fit(line, avg10['XT'], avg10['D47'], sigma=avg10['D47_sterr'])
popt5, pcov5 = curve_fit(line, avg10['XT'], avg10['D47_avg'])
popt, pcov = curve_fit(line, avg['XT'], avg['D47_avg'])
poptt, pcovt = curve_fit(line, travertine['XT'], travertine['D47'])
poptl, pcovl = curve_fit(line, leeds['XT'], leeds['D47'])
poptE, pcovE = curve_fit(line, ETH['XT'], ETH['D47'])
poptc, pcovc = curve_fit(line, cavepearls['XT'], cavepearls['D47'])
poptBF, pcovBF = curve_fit(line, magali['XT'], magali['D47'])
popt18, pcov18 = curve_fit(line, avg10['d18O_avg'], avg10['D47_avg'])
popt18b, pcov18b = curve_fit(line, avg10['T'], avg10['d18O_avg'] )
# poptb, pcovb = curve_fit(line, big['XT'], big['D47_avg'])
# poptnb, pcovnb = curve_fit(line, notbig['XT'], notbig['D47_avg'])
popts, pcovs = curve_fit(line, site['XT'], site['D47_avg'])
str_all = 'y='+ '%.3f' % popt_all[0] + 'x+' + '%.3f' % popt_all[1]
strc = 'y='+ '%.3f' % popt5[0] + 'x+' + '%.3f' % popt5[1]
strt = 'y='+ '%.3f' % poptt[0] + 'x+' + '%.3f' % poptt[1]
strl = 'y='+ '%.3f' % poptl[0] + 'x+' + '%.3f' % poptl[1]
strE = 'y='+ '%.3f' % poptE[0] + 'x+' + '%.3f' % poptE[1]
strcp = 'y='+ '%.3f' % poptc[0] + 'x+' + '%.3f' % poptc[1]
# strb = 'y='+ '%.3f' % poptb[0] + 'x+' + '%.3f' % poptb[1]
# strnb = 'y='+ '%.3f' % poptnb[0] + 'x+' + '%.3f' % poptnb[1]
strs = 'y='+ '%.4f' % popts[0] + 'x+' + '%.3f' % popts[1]

strB = 'y=0.0422x+0.208'

# r squared
residuals = avg10['D47_avg']- line(avg10['XT'], popt5[0], popt5[1])
ss_res = np.sum(residuals**2)
ss_tot = np.sum((avg10['D47_avg']-np.mean(avg10['D47_avg']))**2)
r_squared = 1 - (ss_res / ss_tot)
strr = 'r$^2$=' '%.2f' % r_squared

residuals_reps = reps['D47'] - line(reps['XT'], popt_all[0], popt_all[1])
ss_res_reps = np.sum(residuals_reps**2)
ss_tot_reps = np.sum((reps['D47']-np.mean(reps['D47']))**2)
r_squared_reps = 1 - (ss_res_reps / ss_tot_reps)
strr_reps = 'r$^2$=' '%.2f' % r_squared_reps

residuals_site = site['D47_avg'] - line(site['XT'], popts[0], popts[1])
ss_res_site = np.sum(residuals_site**2)
ss_tot_site = np.sum((site['D47_avg']-np.mean(site['D47_avg']))**2)
r_squared_site = 1 - (ss_res_site / ss_tot_site)
strr_site = 'r$^2$=' '%.2f' % r_squared_site


# confidence interval
xplot = linspace(np.min(avg10['XT'])-0.1, np.max(avg10['XT'])+0.1, 100)
ps = np.random.multivariate_normal(popt5, pcov5, 10000)
ysample = np.asarray([line(xplot, *pi) for pi in ps])
lower = percentile(ysample, 2.5, axis = 0)
upper = percentile(ysample, 97.5, axis = 0)

ps_all = np.random.multivariate_normal(popt_all, pcov_all, 10000)
ysample_all = np.asarray([line(xplot, *pi) for pi in ps_all])
lower_all = percentile(ysample_all, 2.5, axis = 0)
upper_all = percentile(ysample_all, 97.5, axis = 0)

ps_site = np.random.multivariate_normal(popts, pcovs, 10000)
ysample_site = np.asarray([line(xplot, *pi) for pi in ps_site])
lower_site = percentile(ysample_site, 2.5, axis = 0)
upper_site = percentile(ysample_site, 97.5, axis = 0)

ps_bonifacie = np.random.multivariate_normal(poptBF, pcovBF, 10000)
ysample_BF = np.asarray([line(xplot, *pi) for pi in ps_bonifacie])
lower_BF = percentile(ysample_BF, 2.5, axis= 0)
upper_BF = percentile(ysample_BF, 97.5, axis = 0)



# plot my data plus other peoples basically the OG plot
# fig1 = plt.figure()
# plt.rc('font', family='Helvetica')
# ax = fig1.add_subplot(1,1,1)
# ax.plot(reps['XT'], reps['D47'], '+', color = '0.75', label = 'reps',  zorder=1)
# ax.errorbar(avg['XT'], avg['D47_avg'], yerr = avg['D47_sterr'],
#             fmt = 'o',  color = 'b', ecolor= 'b', label = 'avg')
# ax.errorbar(avg10['XT'], avg10['D47_avg'], yerr = avg10['D47_sterr'],
#             fmt = 'o', color = 'm', ecolor = 'm', label = 'avg>10')
# ax.plot(travertine['XT'], travertine['D47'], 'o', color = 'g', label = 'travertine')
# ax.plot(leeds['XT'], leeds['D47'], 'o', color = 'r', label = 'leeds')
# ax.plot(ETH['XT'], ETH['D47'], 'o', color = 'k', label = 'ETH')
# ax.plot(cavepearls['XT'], cavepearls['D47'], 'o', color = 'y', label = 'cavepearls')
# ax.plot(xplot, line(xplot, popt[0], popt[1]), 'k--', linewidth=0.5)
# ax.plot(xplot, line(xplot, popt5[0], popt5[1]), 'm-')
# ax.plot(xplot, line(xplot, poptt[0], poptt[1]), 'g--')
# ax.plot(xplot, line(xplot, poptl[0], poptl[1]), 'r--')
# ax.plot(xplot, line(xplot, poptE[0], poptE[1]), 'k--')
# ax.plot(xplot, line(xplot, poptc[0], poptc[1]), 'y--')
#
# #ax.plot(xplot, lower, 'm--')
# #ax.plot(xplot, upper, 'm--')
# ax.fill_between(xplot, upper, lower, facecolor='m', alpha=0.5)
#
# txt = strc + '\n' + strr + '\n' + strt + '\n' + strl + '\n' + strE + '\n' + strcp
# # txt = strc + '\n' + strr
# ax.text(0.05, 0.72, txt, transform=ax.transAxes)
# ax.set_xlabel('$10^6/$T$^2 ($K$)$')
# ax.set_ylabel('$\Delta_{47}$ (\u2030)')
# ax.legend(loc='lower right')
# ax.set_xlim([11.5, 13.75])
# ax.set_ylim([0.65, 0.85])
#
# # a2 = fig1.add_subplot(2,1,2)
# # a2.plot(avg10['XT'], residuals, 'ko')
# # a2.set_xlabel('$10^6/$T$^2 ($K$)$')
# # a2.set_ylabel('residual')
# #plt.show()
# #avg.to_csv('avgout.csv')
# plt.savefig('12_7_all_noerr.pdf')

# plot just my data
# fig2 = plt.figure()
# plt.rc('font', family='Helvetica')
# ax2 = fig2.add_subplot(1,1,1)
# ax2.plot(reps['XT'], reps['D47'], '+', color = '0.75', label = 'Replicates',  zorder=1)
# ax2.errorbar(avg10['XT'], avg10['D47_avg'], yerr = avg10['D47_sterr'],
#             fmt = 'o', color = 'm', ecolor = 'm', label = 'Sample Averages')
# # plot regression based on Averages
# ax2.plot(xplot, line(xplot, popt5[0], popt5[1]), 'm-')
# ax2.fill_between(xplot, upper, lower, facecolor='m', alpha=0.5)
# txt2 = strc + '\n' + strr
#
# # plot regression based on Replicates
# # ax2.plot(xplot, line(xplot, popt_all[0], popt_all[1]), 'm-')
# # ax2.fill_between(xplot, upper_all, lower_all, facecolor='m', alpha=0.5)
# # txt2 = str_all + '\n' + strr_reps
# ax2.text(0.47, 0.02, txt2, transform=ax.transAxes)
# ax2.set_xlabel('$10^6/$T$^2 ($K$)$')
# ax2.set_ylabel('$\Delta_{47}$ (\u2030)')
# ax2.legend(loc='upper left')
# ax2.set_xlim([11.5, 13.75])
# ax2.set_ylim([0.65, 0.85])
# plt.savefig('12_7_data_noerr.pdf')


# plot regression comparisons
# fig3 = plt.figure()
# plt.rc('font', family='Helvetica')
# ax3 = fig3.add_subplot(1,1,1)
# ax3.errorbar(avg10['XT'], avg10['D47_avg'], yerr = avg10['D47_sterr'],
#             fmt = 'o', color = 'm', ecolor = 'm', label = 'This Study')
# ax3.plot(travertine['XT'], travertine['D47'], 'd', color = 'g', label = 'ETH travertine')
# ax3.plot(leeds['XT'], leeds['D47'], '^', color = 'r', label = 'ETH synthetic')
# # ax3.plot(ETH['XT'], ETH['D47'], 's', color = 'k', label = 'ETH synthetic')
# # ax3.plot(cavepearls['XT'], cavepearls['D47'], 'H', color = 'y', label = 'Cambridge')
# ax3.plot(xplot, line(xplot, popt5[0], popt5[1]), 'm-')
# ax3.plot(xplot, line(xplot, poptt[0], poptt[1]), 'g--')
# ax3.plot(xplot, line(xplot, poptl[0], poptl[1]), 'r--')
# # ax3.plot(xplot, line(xplot, poptE[0], poptE[1]), 'k--')
# # ax3.plot(xplot, line(xplot, poptc[0], poptc[1]), 'y--')
# ax3.plot(xplot, line(xplot, 0.0422, 0.2082),      'c--', label = 'Bonifacie et al 2017')
# ax3.fill_between(xplot, upper, lower, facecolor='m', alpha=0.5)
# txt3 = strt + '\n' + strl + '\n' + strB + '\n' + strc
# ax3.text(-0.01, 0.55, txt3, transform=ax.transAxes)
# ax3.set_xlabel('$10^6/$T$^2 ($K$)$')
# ax3.set_ylabel('$\Delta_{47}$ (\u2030)')
# ax3.legend(loc='lower right')
# ax3.set_xlim([11.5, 13.75])
# ax3.set_ylim([0.65, 0.85])
# plt.savefig('12_7_compare_noerr.pdf')

# # plot regression comparisons
# fig3 = plt.figure()
# plt.rc('font', family='Helvetica')
# ax3 = fig3.add_subplot(1,1,1)
# ax3.errorbar(avg10['XT'], avg10['D47'], yerr = avg10['D47_sterr'],
#             fmt = 'o', color = 'm', ecolor = 'm', label = 'This Study')
# ax3.plot(travertine['XT'], travertine['D47'], 'd', color = 'g', label = 'ETH travertine')
# ax3.plot(leeds['XT'], leeds['D47'], '^', color = 'r', label = 'ETH synthetic')
# # ax3.plot(ETH['XT'], ETH['D47'], 's', color = 'k', label = 'ETH synthetic')
# # ax3.plot(cavepearls['XT'], cavepearls['D47'], 'H', color = 'y', label = 'Cambridge')
# xplot2 = linspace(8-0.1, 14+0.1, 100)
# ax3.plot(xplot, line(xplot, popt5[0], popt5[1]), 'm-')
# ax3.plot(xplot2, line(xplot2, 0.0464, 0.1535), 'g--')
# ax3.plot(xplot2, line(xplot2, 0.0464, 0.132), 'r--')
# # ax3.plot(xplot, line(xplot, poptE[0], poptE[1]), 'k--')
# # ax3.plot(xplot, line(xplot, poptc[0], poptc[1]), 'y--')
# ax3.plot(xplot2, line(xplot2, 0.0422, 0.2082),      'c--', label = 'Bonifacie et al 2017')
# ax3.fill_between(xplot, upper, lower, facecolor='m', alpha=0.5)
# txt3 = strt + '\n' + strl + '\n' + strB + '\n' + strc
# ax3.text(-0.01, 0.55, txt3, transform=ax.transAxes)
# ax3.set_xlabel('$10^6/$T$^2 ($K$)$')
# ax3.set_ylabel('$\Delta_{47}$ (\u2030)')
# ax3.legend(loc='best')
# ax3.set_xlim([8, 14])
# ax3.set_ylim([0.45, 0.75])
# plt.savefig('12_7_compare_noerr.eps', format = 'eps', dpi = 1000)


# plot species specific
fig4 = plt.figure()
plt.rc('font', family='Helvetica')
ax4 = fig4.add_subplot(1,1,1)
ax4b = ax4.twiny()
# ax4.plot(xplot, line(xplot, popt5[0], popt5[1]), 'k--', label = 'This Study')
# ax4.fill_between(xplot, upper, lower, facecolor='k', alpha=0.2)
ax4.plot(xplot, line(xplot, popts[0], popts[1]), 'k--', label = 'This Study' )
ax4.fill_between(xplot, upper_site, lower_site, facecolor='k', alpha=0.2)
ax4.errorbar(cibs['XT'], cibs['D47_avg'], yerr = cibs['D47_sterr'],
            fmt = 'o', color = 'r', ecolor = 'r', label = 'C. pachyderma')
ax4.errorbar(lent['XT'], lent['D47_avg'], yerr = lent['D47_sterr'],
            fmt = '^', color = 'g', ecolor = 'g', label = 'Lenticulina spp')
ax4.errorbar(pyrgo['XT'], pyrgo['D47_avg'], yerr = pyrgo['D47_sterr'],
            fmt = 'h', color = 'c', ecolor = 'c', label = 'Pyrgo spp')
ax4.errorbar(uvi['XT'], uvi['D47_avg'], yerr = uvi['D47_sterr'],
            fmt = 's', color = 'b', ecolor = 'b', label = 'U. mediteranea')
ax4.errorbar(elegans['XT'], elegans['D47_avg'], yerr = elegans['D47_sterr'],
            fmt = 'p', color = 'darkviolet', ecolor = 'darkviolet', label = 'H. elegans')
ax4.errorbar(other['XT'], other['D47_avg'], yerr = other['D47_sterr'],
            fmt = 'd', color = 'm', ecolor = 'm', label = 'Assorted')
# ax4.plot(xplot, line(xplot, 0.0422, 0.2082),      'c--', label = 'Bonifacie et al 2017')
ax4.set_xlabel('$10^6/$T$^2 ($K$)$')
ax4.set_ylabel('$\Delta_{47}$ (\u2030)')
ax4.legend(loc='lower right')
ax4.set_xlim([11.55, 13.7])
ax4.set_ylim([0.65, 0.85])
ax4Ticks = ax4.get_xticks()
ax4bTicks = ax4Ticks
def tick_function(X):
    V = -273.15 + ((1e6/X)**(0.5))
    return ['%.1f' %z for z in V]
ax4b.set_xticks(ax4bTicks)
ax4b.set_xbound(ax4.get_xbound())
ax4b.set_xticklabels(tick_function(ax4bTicks))
ax4b.set_xlabel('T ($^\circ$C)')
plt.savefig('02_13_18_species_noerr.pdf')

# d18O plotz

fig5 = plt.figure()
plt.rc('font', family= 'Helvetica')
ax5= fig5.add_subplot(1,1,1)
x18plot = linspace(np.min(avg10['d18O_avg'])-0.1, np.max(avg10['d18O_avg'])+0.1, 100)
# ax5.errorbar(avg10['d18O_avg'], avg10['D47_avg'], yerr = avg10['D47_sterr'], xerr = avg10['d18O_sterr'],
            # fmt = 'd', color = 'k', ecolor = 'k', label = 'Averages >10')
ax5.errorbar(cibs['d18O_avg'], cibs['D47_avg'], yerr = cibs['D47_sterr'], xerr = cibs['d18O_sterr'],
            fmt = 'o', color = 'r', ecolor = 'r', label = 'C. pachyderma')
ax5.errorbar(lent['d18O_avg'], lent['D47_avg'], yerr = lent['D47_sterr'], xerr = lent['d18O_sterr'],
            fmt = '^', color = 'g', ecolor = 'g', label = 'Lenticulina spp')
ax5.errorbar(pyrgo['d18O_avg'], pyrgo['D47_avg'], yerr = pyrgo['D47_sterr'], xerr = pyrgo['d18O_sterr'],
            fmt = 'h', color = 'c', ecolor = 'c', label = 'Pyrgo spp')
ax5.errorbar(uvi['d18O_avg'], uvi['D47_avg'], yerr = uvi['D47_sterr'], xerr = uvi['d18O_sterr'],
            fmt = 's', color = 'b', ecolor = 'b', label = 'U. mediteranea')
ax5.errorbar(elegans['d18O_avg'], elegans['D47_avg'], yerr = elegans['D47_sterr'], xerr = elegans['d18O_sterr'],
            fmt = 'p', color = 'darkviolet', ecolor = 'darkviolet', label = 'H. elegans')
ax5.errorbar(other['d18O_avg'], other['D47_avg'], yerr = other['D47_sterr'], xerr = other['d18O_sterr'],
            fmt = 'd', color = 'm', ecolor = 'm', label = 'Assorted')
ax5.plot(x18plot, line(x18plot, popt18[0], popt18[1]), 'k--')
ax5.legend(loc = 'upper left')
ax5.set_xlabel('$\delta^{18}$O (\u2030)')
ax5.set_ylabel('$\Delta_{47}$ (\u2030)')
plt.savefig('02_13_18_d18O.pdf')

# site average plot
fig6 = plt.figure()
plt.rc('font', family = 'Helvetica')
ax6 = fig6.add_subplot(1,1,1)
ax6b = ax6.twiny()
ax6.plot(xplot, line(xplot, popts[0], popts[1]), 'k--', label = 'This Study' )
ax6.fill_between(xplot, upper_site, lower_site, facecolor='k', alpha=0.2)
# ax6.plot(xplot, line(xplot, 0.0422, 0.2082),      'c--', label = 'Bonifacie et al 2017')
ax6.errorbar(site['XT'], site['D47_avg'], yerr = site['D47_sterr'],
            fmt = 'o', color = 'k', ecolor = 'k', label = 'Site Averages', zorder=10)
ax6.plot(reps['XT'], reps['D47'], '+', color = '0.75', label = 'Replicates',  zorder=1)
ax6.plot(xplot, line(xplot, 0.0422, 0.2082),      'c--', label = 'Bonifacie et al 2017')
ax6.plot(magali['XT'], magali['D47'], 'd', color = 'c', alpha=0.2, label = 'Bonifacie et al 2017')

ax6.set_xlabel('$10^6/$T$^2 ($K$)$')
ax6.set_ylabel('$\Delta_{47}$ (\u2030)')
ax6.legend(loc='upper left')
ax6.set_xlim([11.55, 13.7])
ax6.set_ylim([0.65, 0.85])
txts = strs + '\n' + strr_site
ax6.text(0.77, 0.02, txts, transform=ax6.transAxes)
ax6Ticks = ax6.get_xticks()
ax6bTicks = ax6Ticks
ax6b.set_xticks(ax6bTicks)
ax6b.set_xbound(ax6.get_xbound())
ax6b.set_xticklabels(tick_function(ax6bTicks))
ax6b.set_xlabel('T ($^\circ$C)')
# ax6b.set_xlim([-1.5, 20])
plt.savefig('02_13_18_site_noerr.pdf')

# different d18O plot
fig7 = plt.figure()
plt.rc('font', family = 'Helvetica')
ax7 = fig7.add_subplot(1,1,1)
x18plot2 = linspace(np.min(avg10['T'])-1, np.max(avg10['T'])+1, 100)
# ax5.errorbar(avg10['d18O_avg'], avg10['D47_avg'], yerr = avg10['D47_sterr'], xerr = avg10['d18O_sterr'],
            # fmt = 'd', color = 'k', ecolor = 'k', label = 'Averages >10')
ax7.errorbar(cibs['T'], cibs['d18O_avg'], yerr = cibs['d18O_sterr'],
            fmt = 'o', color = 'r', ecolor = 'r', label = 'C. pachyderma')
ax7.errorbar(lent['T'], lent['d18O_avg'], yerr = lent['d18O_sterr'],
            fmt = '^', color = 'g', ecolor = 'g', label = 'Lenticulina spp')
ax7.errorbar(pyrgo['T'], pyrgo['d18O_avg'], yerr = pyrgo['d18O_sterr'],
            fmt = 'h', color = 'c', ecolor = 'c', label = 'Pyrgo spp')
ax7.errorbar(uvi['T'], uvi['d18O_avg'], yerr = uvi['d18O_sterr'],
            fmt = 's', color = 'b', ecolor = 'b', label = 'U. mediteranea')
ax7.errorbar(elegans['T'], elegans['d18O_avg'], yerr = elegans['d18O_sterr'],
            fmt = 'p', color = 'darkviolet', ecolor = 'darkviolet', label = 'H. elegans')
ax7.errorbar(other['T'], other['d18O_avg'], yerr = other['d18O_sterr'],
            fmt = 'd', color = 'm', ecolor = 'm', label = 'Assorted')
ax7.plot(x18plot2, line(x18plot2, popt18b[0], popt18b[1]), 'k--')
ax7.legend(loc = 'lower left')
ax7.set_xlabel('T ($^\circ$C)')
ax7.set_ylabel('$\delta^{18}$O (\u2030)')
plt.savefig('02_13_18_d18OvT.pdf')

# DATA COMPARE
fig8 = plt.figure()
plt.rc('font', family = 'Helvetica')
ax8 = fig8.add_subplot(1,1,1)
ax8b = ax8.twiny()
ax8.plot(xplot, line(xplot, popts[0], popts[1]), 'k--', label = 'This Study' )
ax8.fill_between(xplot, upper_site, lower_site, facecolor='k', alpha=0.2)
ax8.plot(xplot, line(xplot, 0.0422, 0.2082),      'c--', label = 'Bonifacie et al 2017')
ax8.errorbar(site['XT'], site['D47_avg'], yerr = site['D47_sterr'],
            fmt = 'o', color = 'k', ecolor = 'k', label = 'Site Averages')
# ax8.plot(reps['XT'], reps['D47'], '+', color = '0.75', label = 'Replicates',  zorder=1)
# ax8.plot(magali['XT'], magali['D47'], 'o', color = 'blue', label = 'Bonifacie et al 2017')
ax8.fill_between(xplot, upper_BF, lower_BF, facecolor='c', alpha=0.2)

ax8.set_xlabel('$10^6/$T$^2 ($K$)$')
ax8.set_ylabel('$\Delta_{47}$ (\u2030)')
ax8.legend(loc='upper left')
# ax8.set_xlim([11.55, 13.7])
# ax8.set_ylim([0.65, 0.85])
txts = strs + '\n' + strr_site
ax8.text(0.77, 0.02, txts, transform=ax8.transAxes)
ax8Ticks = ax8.get_xticks()
ax8bTicks = ax8Ticks
ax8b.set_xticks(ax8bTicks)
ax8b.set_xbound(ax8.get_xbound())
ax8b.set_xticklabels(tick_function(ax8bTicks))
ax8b.set_xlabel('T ($^\circ$C)')
# ax8b.set_xlim([-1.5, 20])
plt.savefig('02_13_18_compare_magali_noerr.pdf')
