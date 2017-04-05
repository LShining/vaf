from hoomd_script import *
import numpy as np
import scipy
import math
import pickle
import matplotlib
from datetime import datetime
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal
from vafnewpath import *


context.initialize()

def calTau(length,epsilon,sigma):
	test = vcf()
	test.init(periodLength = 1400, system = system, printFigure = 1, numOfParticles = 300)
 
#run 1,000,000 time steps

	run(5000)
	run(length, callback_period=10,callback=test.main)

#process and output the figure
	return test.dataProcess()


def convergeTau(length = 16000, epsilon=1.0, sigma=1.0):
	if epsilon != 1.0:
		lj.pair_coeff.set('A', 'A', epsilon=epsilon)
	if sigma != 1.0:
		lj.pair_coeff.set('A', 'A', sigma=sigma)
	length = length
	lastTau = calTau(length , epsilon , sigma)
	while True:
		length += 1000
		tau = calTau(length, epsilon , sigma)
		if abs(tau-lastTau)/lastTau < 0.005:
			return tau
		lastTau = tau


# create random particles of name A
system = init.create_random(N=1000, phi_p=0.01, name='A', seed=4864)
# specify Lennard-Jones interactions between particle pairs
lj = pair.lj(r_cut=3.0)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
dump.dcd(filename=str(datetime.now())+'trajectory.dcd', period=10)
update.zero_momentum()

# integrate at constant temperature
all = group.all();
integrate.mode_standard(dt=0.0025)
integrate.nvt(group=all, T=1.4,tau=0.5)


epsilonList = np.arange(0.6,1.4,0.1)
tauList = [ convergeTau(15000,x,1.0) for x in epsilonList]

font = {'family': 'serif',  'color':  'black',   'weight': 'normal',   'size': 12 ,  }
plt.title('Tau against epsilon', fontdict=font)
plt.xlabel('epsilon', fontdict=font)
plt.ylabel('tau', fontdict=font)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(epsilonList,tauList)
fig.savefig('TauAndEpsilon.png')
print("Calculation complete")

