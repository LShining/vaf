from hoomd_script import *
from datetime import datetime
import numpy as np
import scipy
import matplotlib
from scipy.optimize import curve_fit
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal
import collections
from numba import autojit
class vcf:
# poi?
#------------variables------------------------------
	vaf = []			#This is where the final velocity autocorrelation functions will go
	velocitySeries = []           	#This is to store the particle velocities of all time
	productSeries = []            	#This is to store the products of all particle velocities
	tau = 0 			#Decay time
	nperiod = 0             	#Number of the subperiods calculated
	periodLength = 0		#Value of largest time interval used
	numOfParticles = 0		#Number of all particles in the simulation
	correlationDemension = 0	#'full' : calculate over all demensions; 'x' : calculate x-direction properties only
	correlationType = 0		#Can be 'velocity', 'orientation', or 'norm'
	callback_period = 0		#How often is a subperoid started, 0 to start another immediately after one is finished
	printFigure = 0			#Flag | 1: print a figure of the curve
	system = 0			#Pointer for system in hoomd to be passed
	velocityRecorder = 0		#Dummy function pointer that takes effect after initialization
	correlationCalculator = 0	#Same as above
	tau = 0				#Decay time in unit of steps
#------------functions-----------------------------

	def norm(self, somelist):
		return sum([x*x for x in somelist])**0.5
	def normalizeVector(self, somelist):
		return [x/self.norm(somelist) for x in somelist]
	
	#Recorder functions called by vcf_main
	
	def vRecorderFull(self):
		for particleLabel in range(0,self.numOfParticles):
			self.velocitySeries[particleLabel].append(self.system.particles[particleLabel].velocity)
	def vRecorderX(self):
		for particleLabel in range(0,self.numOfParticles):
			self.velocitySeries[particleLabel].append(system.particles[particleLabel].velocity[0])
	def nRecorderFull(self):
		for particleLabel in range(0,self.numOfParticles):
			self.velocitySeries[particleLabel].append(self.norm(system.particles[particleLabel].velocity))
	def nRecorderX(self):
		for particleLabel in range(0,self.numOfParticles):
			self.velocitySeries[particleLabel].append(abs(system.particles[particleLabel].velocity[0]))
	def oRecorderFull(self):
		for particleLabel in range(0,self.numOfParticles):
			self.velocitySeries[particleLabel].append(self.normalizeVector(system.particles[particleLabel].velocity))
	def oRecorderX(self):
		for particleLabel in range(0,self.numOfParticles):
			self.velocitySeries[particleLabel].append(np.sign(system.particles[particleLabel].velocity[0]))
	
	#Calulation functions called by vcf_main
	
	def scalarCorrelationCalculator(self):
		for particleVelocitySerie in self.velocitySeries:
			corrp = [item*particleVelocitySerie[0] for item in particleVelocitySerie]
			self.productSeries.append(corrp)   #append array to pproduct as a list
		for periodIndex in range(0 , self.periodLength):   #walk through all time intervals
			sumDummy = 0
			for particleLabel in range(0,self.numOfParticles):    #walk through all particles
				sumDummy += self.productSeries[particleLabel][periodIndex]
			sumDummy /= self.numOfParticles   #take average over all particles
			try:
				self.vaf[periodIndex] += sumDummy             #the average is stored into vaf
			except IndexError:			#break when exceeding the limit
				break
	@autojit
	def vectorCorrelationCalculator(self):
		for vectorDemension in range(0,3): #calculate 3 demensions individually then sum up
			for particleVelocitySerie in self.velocitySeries:
				singleDemensionVelocitySerie = [ everyParticleVelocity[vectorDemension] for everyParticleVelocity in particleVelocitySerie]
#				corrp = [item*singleDemensionVelocitySerie[0] for item in singleDemensionVelocitySerie]
				corrp = [1 for item in singleDemensionVelocitySerie]
				self.productSeries.append(corrp)
			for periodIndex in range(0,self.periodLength): 
				sumDummy = 0
				for particleLabel in range(0,self.numOfParticles):
					try:
						sumDummy += self.productSeries[particleLabel][periodIndex]
					except IndexError:
						print("Label",particleLabel,"Index",periodIndex)
						break
				sumDummy /= self.numOfParticles
				try:
					self.vaf[periodIndex] += sumDummy
				except IndexError:	
					break	
			self.productSeries = []
	#Initialization function, supposed to be called before 'run'
	
	def __init__(self, system, periodLength = 1000, numOfParticles = 1000,  correlationDemension = 'full', correlationType = 'velocity',  printFigure = 1 , callback_period= 10):
		self.periodLength = periodLength
		self.numOfParticles = numOfParticles
		self.correlationDemension = correlationDemension
		self.correlationType = correlationType
		self.vaf = [0] * periodLength
		self.velocitySeries = [ collections.deque(maxlen = periodLength)  for x in range(0,self.numOfParticles) ] 
		self.system = system
		self.printFigure = printFigure
		self.callback_period = callback_period
		if self.correlationType == 'velocity' :
			if self.correlationDemension == 'full':
				self.velocityRecorder = self.vRecorderFull
				self.correlationCalculator = self.vectorCorrelationCalculator
			elif self.correlationDemension == 'x':
				self.velocityRecorder = self.vRecorderX
				self.correlationCalculator = self.scalarCorrelationCalculator
			else:
				raise TypeError("Can only accept 'full' or 'x' as correlationDemension so far")
		elif self.correlationType == 'norm' :
			if self.correlationDemension == 'full':
				self.velocityRecorder = self.periodLengthnRecorderFull
				self.correlationCalculator = self.scalarCorrelationCalculator
			elif self.correlationDemension == 'x':
				self.velocityRecorder = self.nRecorderX
				self.correlationCalculator = self.scalarCorrelationCalculator
			else:
				raise TypeError("Can only accept 'full' or 'x' as correlationDemension so far")
		elif self.correlationType == 'orientation':
			if self.correlationDemension == 'full':
				self.velocityRecorder = self.oRecorderFull
				self.correlationCalculator = self.vectorCorrelationCalculator
			elif self.correlationDemension == 'x':
				self.velocityRecorder = self.oRecorderX
				self.correlationCalculator = self.scalarCorrelationCalculator
			else:
				raise TypeError("Can only accept 'full' or 'x' as correlationDemension so far")
		else:
			raise TypeError("Can only accept 'velocity', 'norm', or 'orientation' as correlationDemension so far")

#Main function; records the real-time x-velocities of all N particles and calculate VAF for every subperoid

	def main(self , cur_step):
		self.velocityRecorder()
		if len(self.velocitySeries[0]) >= self.periodLength:
			self.nperiod += 1
			self.correlationCalculator()
	def dataProcess(self):
		self.vaf = [x / self.vaf[0] for x in self.vaf] #normalization
		regressionScope = 0.5*self.periodLength
		for value in self.vaf:
			if value < 0.2:
				regressionScope = self.vaf.index(value)
				break
		regressionScope = int(regressionScope)
		vdata = self.vaf[:regressionScope]
		tdata = range(1,regressionScope+1)
		vdata = np.array(vdata)
		tdata = np.array(tdata)
		def func(t, tau):
			return np.exp(-t / tau)
		self.tau,_ = curve_fit(func, tdata, vdata)
		print("current tau:",self.tau)
		if self.printFigure == 1:
			self.printFig()
		return self.tau
	def printFig(self):
		font = {'family': 'serif',  'color':  'black',   'weight': 'normal',   'size': 12 ,  }
		xaxis = range(0, self.periodLength)
		plt.title('Velocity Autocorrelation Function', fontdict=font)
		plt.xlabel('time(steps)', fontdict=font)
		plt.ylabel('Normalized VAcF', fontdict=font)
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(xaxis , self.vaf)
		fig.savefig(str(datetime.now())+'.png')
	def systemDensity(self):
		return self.numOfParticles/((86.8233 * 2) ** 3)

