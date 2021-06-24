from neuron import h
from neuron.units import ms, mV
import numpy as np


def setupsoma(modelparams, printout=False):


	# Create the compartment and insert passive
	soma = h.Section(name='soma')
	# Passive Channels
	soma.insert('pas')

	# Set the passive and leak parameters
	soma.cm = modelparams.cm # uF/cm2
	soma.g_pas = modelparams.gl* 1e-3 # mS/cm2  -> NEURON model: S/cm2 
	soma.e_pas = modelparams.ErL # mV
        
	# Load the channels
	if modelparams.has_vchan: 
		h.load_file(modelparams.channelsesfiles)
	else:
		if printout: # Check model
			print("\t**Passive Model**")
			print("\t-Resting potential and init E0 {:.4f} mV".format(modelparams.E0))
			print("\t-Leak reversal potential {:.4f} mV".format(soma.e_pas))
			print("\t-Leak conductance {:.4f} mS/cm2".format(soma.g_pas* 1e3))
			print("\t-Capa {:.4f} uF/cm2".format(soma.cm))
			
		return soma, None, None, None

	# Geometry
	if modelparams.geo: # we are stochastic 
		soma.diam = modelparams.geo['diam'] # um   
		soma.L = modelparams.geo['L']

		# Insert the voltage-dep conductances
		# we are stochastic : Insert the voltage-dep conductances as pt process

		# Sodium (custom)
		stch_na = h.stch_ctr_na(soma(0.5))
		[stch_na.gmax, stch_na.Nsingle] = modelparams.get_gmax_numch_stoch('na',soma(0.5).area())
		# Delayed rectifier Potassium (custom)
		stch_kdr = h.stch_ctr_kdr(soma(0.5))
		[stch_kdr.gmax, stch_kdr.Nsingle] = modelparams.get_gmax_numch_stoch('kdr',soma(0.5).area())
        
		# Additional sub-threshold channel (NOISY+) (custom)
		stch_subchan = h.stch_ctr_subchan(soma(0.5))
		[stch_subchan.gmax, stch_subchan.Nsingle] = modelparams.get_gmax_numch_stoch('subchan',soma(0.5).area())
                
	else: # full deterministic model 
		# Sodium (custom)
		soma.insert('ch_ctr_na')
		# Delayed rectifier Potassium (custom)
		soma.insert('ch_ctr_kdr')            
		# Additional sub-threshold channel (NOISY+) (custom)
		soma.insert('ch_ctr_subchan')       
		# Modify Channels max conductances
		list_ch = ['na','kdr','subchan']
		for ch in list_ch:
			u=soma(0.5).__getattribute__('ch_ctr_'+ch)
			u.gmax = modelparams.gx[ch] * 1e-3 # mS/cm2  -> NEURON model: S/cm2
		
		stch_na=None
		stch_kdr=None
		stch_subchan=None

	# Setup the v-dep channels parameters
	# Reversal potentials - define after insert the point process, otherwise ena and ek won't be recognised
	soma.ena = modelparams.ena # mV
	soma.ek = modelparams.ek # mV

	if printout: # Check model
		print("\t-Channels file: "+modelparams.channelsesfiles)
		print("\t-Resting potential and init E0 {:.4f} mV".format(modelparams.E0))
		print("\t-Leak reversal potential {:.4f} mV".format(soma.e_pas))
		print("\t-Leak conductance {:.4f} mS/cm2".format(soma.g_pas* 1e3))
		print("\t-Capa {:.4f} uF/cm2".format(soma.cm))
        
		print("\t-AP threshold {:.4f} mV".format(modelparams.APthreshold))
        
		if modelparams.geo:
			print("\t-Surface Area {:.4f} um2".format(soma(0.5).area()))
			print(f'\t-Number of Na channels {stch_na.Nsingle} (0 means run deterministic)')
			print(f'\t-Number of Kdr channels {stch_kdr.Nsingle} (0 means run deterministic)')
			print(f'\t-Number of SubT channels {stch_subchan.Nsingle} (0 means run deterministic)')
			#print("\t-Number of Na channels "+str(stch_na.Nsingle)+)
			#print("\t-Number of Kdr channels "+str(stch_kdr.Nsingle))

	#print(ch+" gbar {:.4f} mS/cm2".format(soma(0.5).__getattribute__(channelrepr[key_group]+ch).gmax* 1e3))
	#print(h.units('g_'+channelrepr[key_group]+ch))

	return soma, stch_na, stch_kdr, stch_subchan


class Synapse:

	def __init__(self, soma, simparams):

		self.intype = simparams['intype']
		if (simparams['intype']=='isyn') or (simparams['intype']=='iclamp'):
			if simparams['syntype']=='isyn' :
				t_i = np.linspace(0.0, simparams['simulation_length'], int(simparams['simulation_length']/h.dt))
				self._setup_isyn(soma, simparams['synparams'], t_i)
			else:
				# simple IClamp
				self._setup_iclamp(soma, simparams['stimdelay'], simparams['stimduration'])

			self.unit_convert = 1e3*soma(0.5).area()*1e-8 # nA  # SA in um2 convert to cm2
            # uA/cm2 -> nA/cm2 -> nA # SA in um2 convert to 
		elif simparams['intype']=='gsyn':
			self._setup_gsyn(soma, simparams['synparams'])
            # Set the conductance conversion 
			self.unit_convert = soma(0.5).area()*1e-8 # uS  # SA in um2 convert to cm2

		# Check if provided noise
		if 'syn_noise_sd' in simparams:
			self.syn_noise_sd = simparams['syn_noise_sd']
		else:
			self.syn_noise_sd = 0


	def _setup_gsyn(self, soma, synparams):
		# Insert a synapse with a spikegen
		self.synE = h.Exp2Syn( soma(0.5)  ) # bi-exponential synapse (rise and decay time)
		self.synE.tau1=synparams['tau_rise'] # rise time constant (milliseconds)
		self.synE.tau2=synparams['tau_decay'] # decay time constant (milliseconds)
		self.synE.e = synparams['Er']  # reversal potential
						
		# Provide stimulus to the synapse
		self.stimE = h.NetStim() # Make a new stimulator 
		self.stimE.start = synparams['stimdelay'] # start time of stimulus
		#stimE.interval = 50 # interval between the inputs
		self.stimE.number=1 # number of inputs provided
		self.stimE.noise=0 # how noisy the inputs are (between 0 and 1)

		# Attach it to a synapse 
		# netcon: source, target,threshold (not used here),delay, weight
		self.ncE = h.NetCon(self.stimE,self.synE,0,synparams['delay'],0) 
		#ncE.delay = 1 * ms # can set it later also

	def _setup_isyn(self, soma, synparams, t_i):

		# Simulate synaptic input current waveform
		
		g_i = np.exp(-t_i/synparams['tau_decay'])-np.exp(-t_i/synparams['tau_rise'])

		self.g_i_norm = g_i/g_i.max() 
		self.tvec = h.Vector(t_i)
		
		self.synE_iclamp = h.IClamp(soma(0.5))
		# Next lines necessary to allow play
		self.synE_iclamp.delay = 0.0
		self.synE_iclamp.dur = 1e9

	def _setup_iclamp(self, soma, stimdelay, stimduration):

		self.iclamp = h.IClamp(soma(0.5))
		# Next lines necessary to allow play
		self.iclamp.delay = stimdelay
		self.iclamp.dur = stimduration

	def set_ampmax(self, ampmax):

		# Check if noise
		if self.syn_noise_sd>0:
			ampmax = np.random.normal(ampmax, ampmax*self.syn_noise_sd)
		# (Numpy) Draw random samples from a normal (Gaussian) distribution.
		# This artificially creates trial-to-trial variability
		# NOTE: we multiply stdev by ampmax to scale the noise.
		# .... so make sure the noise is normalised when save it!

		ampmax = ampmax*self.unit_convert
		if self.intype == 'gsyn':
			self.ncE.weight[0] = ampmax

		elif self.intype == 'isyn':

			self.gvec = h.Vector(ampmax * self.g_i_norm) 
			# or : 
			# gvec = h.Vector() 
			# gvec.from_python(g_i)
			# Use play to inset the custom waveform
			self.gvec.play(self.synE_iclamp._ref_amp, self.tvec, True)

		else:
			self.iclamp.amp = ampmax
