import os
import glob
import re
import numpy
import matplotlib
import matplotlib.pyplot as plt
from pylab import *

filelist = glob.glob("/home/nicolas/Documents/logs/2013/06/N2PI/12-6-13/*.txt")
#filelist = glob.glob("/home/nicolas/Documents/logs/2013/06/N2PI/12-6-13/n2_scan18.txt")
print filelist

#fig = plt.figure()
for datafile in filelist:
    source = re.search('([^/]+).txt$', datafile).group(1)
    fig = plt.figure()
    print 'Reading "', datafile, '"'
    f = file(datafile)
    lamda, counts, power, norm = numpy.loadtxt(f, unpack  = True, skiprows = 1)
    data = zip(lamda,counts,power,norm)
    data = [i for i in data if i[2] > 10E-6 and i[0] < 2*237.15]
    lamda = [i[0]/2 for i in data]
    power = [i[2]*1E6 for i in data]
    norm = [i[3]*1E-6 for i in data]
    subplot(111)
    plot (lamda, norm)
    #subplot(212)
    #plot(lamda, power) 
    ax1 = subplot(111)
    lines = {}
    lines['$S_{4/5}$'] = 236.99
    lines['$S_{2/7}$'] = 237.0
    lines['$S_{1/8}$'] = 237.009
    lines['$S_0$'] = 237.022
    lines['$R_{0/4}$'] = (237+0.038)
    lines['$R_5/Q_1$'] = (237+0.049)
    lines['$R_6/Q_2$'] = (237+0.056)
    lines['$P_2$'] = (237+0.072)
    lines['$Q_4$'] = (237+0.078)
    lines['$P_3$'] = (237+0.089)
    lines['$Q_5$'] = (237+0.094)
    lines['$O_3$'] = (237+0.104)
    lines['$P_4$'] = (237+0.110)
    lines['$O_4$'] = 237.135
    dy = -0.4
    P = 0
    for i in sorted(lines.keys(), key = lambda x: lines[x] ):
        if lines[i] > min(lamda) and lines[i] < max(lamda):
            print i
            axvline(lines[i],linewidth = 0.5, ls = '--')
            ax1.annotate(i, (lines[i], 0.96*max(norm)+(1-(-1)**P)*dy), (lines[i] + 0.004, 1*max(norm)+(1-(-1)**P)*dy),
                         arrowprops=dict(color='grey', headwidth = 2, width=0.5, shrink = 0))
            P += 1
    xlabel('$\lambda$ $/nm$')
    ylabel('$counts/pow$')
    grid(color='grey')
#ax2 = subplot(212)
#xlabel('$\lambda$ $/nm$')
#ylabel('$power$')
#grid()
    plt_ = gcf() 
    plt_.set_figwidth(6)
    plt_.set_figheight(3)
    plt_.show()    
    plt_.savefig(source, format = 'png', dpi=180)
    #close()
