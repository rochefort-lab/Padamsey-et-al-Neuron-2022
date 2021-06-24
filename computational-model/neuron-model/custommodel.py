import math

class ModelParameters:

	def __init__(self, key_group, isstochastic, gldrop, APthreshold, chanfile=None, iong=None, clamprest=None, geo=None, ionrunstch=None):

		# Define default values

		# Rest and leak reversal potentials - Units: mV (pre-computed) 
		param_values = {'ctr':{'E0':-75,'ErL':-75},
						'fd':{'E0':-70,'ErL':-70},
						}

		# Default gl
		gl_ctr = 0.1 # mS/cm2 

		self.APthreshold = APthreshold 

		self.cm = 1 # uF/cm2
		self.gl = gldrop*gl_ctr # mS/cm2 
		self.ErL = param_values[key_group]['ErL'] # mV (pre-computed given separation potassium / unspecified leak conductance)
		self.E0 = param_values[key_group]['E0'] # mV (pre-computed) 
        
		self.has_vchan = chanfile!=None and iong!=None # has voltage-dep channels or not

		if self.has_vchan:
			#list_ions = list(iong.keys())
			list_ions = ['na','kdr','subchan']

			if isstochastic:
				self.channelsesfiles = 'stch_'+chanfile+'.ses'
				self.gsingle = {}
				for ion in list_ions:
					self.gsingle[ion] = 20 # pS # Single channel conductance
				if ionrunstch==None:
					self.ionrunstch = {'na':False, 'kdr':False, 'subchan':True}# by default run stochastic only subT chan
				else:
					self.ionrunstch = ionrunstch.copy()
			else:
				self.channelsesfiles = 'ch_'+chanfile+'.ses'
		
			# Reversal potentials 
			self.ena = 55 # mV
			self.ek = -90 # mV
			# active channels ions conductance
			self.gx =  iong.copy() # mS/cm2

		# Save given geometry
		self.geo = geo # um => will be None for deterministic

		# Check if clamp the resting membrane potential; prep a current injection
		self._setup_clamprest(clamprest, param_values)


	def _setup_clamprest(self, clamprest, param_values):
		if clamprest:
			# Compute the current required \\ we simplify by taking the leak conductance only
			if isinstance(clamprest,  str): 
				E0 = param_values[clamprest]['E0']
			else:
				E0 = clamprest
			delta_E = E0 - self.ErL
			self.E0  = E0
			self.iclamprest = self.gl*delta_E # uA/cm2
		else:
			self.iclamprest = None

	def compute_sa(self):
		return math.pi*self.geo*self.geo # um2 
		# / SA (cylinder) = pi*diameter*length = 2*pi*r*h / for us d=h
		# / here = ball \ surface is 4*pi*r^2 , here r=d/2=h/2

	def get_gmax_numch_stoch(self, chan, sa=None):
		"""Compute the number of channels and return max conductance"""
		if sa==None: sa = self.compute_sa()

		gmax0 = self.gx[chan] * 1e3 * 1e-8 # mS/cm2  -> stochastic model uS/um2
		if self.ionrunstch[chan]: # run channel as stochastic
			# Set to Single channel conductance
			gmax = self.gsingle[chan] *1e-6 # pS -> stochastic model: uS            
			N = int( round( (gmax0 / gmax) * sa ) )
		else:  # run channel as deterministic
			gmax = gmax0 * sa # uS 
			N = 0
			# If Nsingle == 0, gating will be deterministic, in which case the parameter gmax will be the maximum conductance of the Point Process when it is completely activated. 
			# https://www.neuron.yale.edu/neuron/static/docs/chanlbild/stochastic/teststoch.html
		return gmax, N

	def compute_currentclamprest(self, sa=None):
		"""For artificial depo/hyperpo"""
		if sa==None: sa = self.compute_sa()

		# Set the current conversion 
		i_convert = 1e3*sa*1e-8 # nA  # SA in um2 convert to cm2
		amp = self.iclamprest*i_convert
		return amp

	
        
def get_simginputparameters(stimdelay, stimduration, intype, synparams=None):

	simginputparams = {}
	simginputparams['intype'] = intype
	simginputparams['stimdelay'] = stimdelay
	simginputparams['stimduration'] = stimduration # stimulation duration
	# Duration of the simulation
	simginputparams['simulation_length'] = stimdelay + stimduration # no additional end-delay because we deal with exp stim
    
	if synparams:
		# Setup synaptic input (only excitatory)
		simginputparams['synparams'] = synparams.copy()
	else:
		# Setup current clamp
		simginputparams['synparams'] = None

	return simginputparams