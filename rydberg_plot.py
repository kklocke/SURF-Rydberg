import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

t = []
m = []
s = []
with open("rydberg_data_lpt03_deltapt2.txt", 'r') as f:
    for line in f:
        line = line.strip().split(' ')
        t.append(float(line[0]))
        m.append(float(line[1]))
        s.append(float(line[2]))

mc = []
with open("rydberg_config_lpt03_deltapt2.txt", 'r') as f:
    for line in f:
        mc.append(line.strip())

titleStr = 'dt: ' + mc[0] + r', $\Delta$: ' + mc[1] + r', $\mu$: ' + mc[2]
titleStr += r', $\lambda$: ' + mc[3] + r', $\gamma$: ' + mc[4] + r', $\rho_0$: ' + mc[5]


plt.errorbar(t, m, yerr=s, elinewidth=0.01, capsize=2., ecolor='red')
plt.yscale('log')
plt.xlabel('T')
plt.ylabel(r'$\langle \rho \rangle$')
plt.title(titleStr)
plt.savefig('rydberg_data_lpt03_deltapt2.png')

def fit_func(x,a,c,d):
    return a*(x**c) + d

def fit_func2(x, m, b):
    return m*x+b

t_min = 150
t_max = 220

xdata = []
ydata = []

for i, tElem in enumerate(t):
    if tElem > t_min and tElem < t_max:
        if m[i] == 0:
            break
        xdata.append(tElem)
        ydata.append(m[i])

xdata = np.array(xdata)
ydata = np.array(ydata)

my_bounds = ([1e19, -14., 0.], [np.inf, -4., 1e-30])
my_bounds2 = ([-20, 0.], [0., np.inf])
plt.figure()
# popt, pcov = curve_fit(fit_func, xdata, ydata, bounds=my_bounds, check_finite=False)
popt, pcov = curve_fit(fit_func2, np.log(xdata), np.log(ydata),bounds=my_bounds2)
# plt.plot(np.log(xdata), np.log(ydata))
# plt.plot(np.log(xdata), fit_func2(np.log(xdata), *popt))
plt.semilogy(xdata, ydata, label='data mean')
plt.semilogy(xdata, np.exp(fit_func2(np.log(xdata), *popt)), 'g--', label="fit")
# plt.semilogy(xdata, fit_func(xdata, *popt), 'g--', label='fit')
# plt.semilogy(xdata, ydata, label='data mean')
plt.xlabel('T')
plt.ylabel(r'$\langle \rho \rangle$')
plt.text(160, 1e-4, r'$\theta = $' + str(popt[0]))
plt.legend()
plt.title(titleStr)
plt.savefig('rydberg_fit_lpt03_deltapt2.png')

print(popt)
