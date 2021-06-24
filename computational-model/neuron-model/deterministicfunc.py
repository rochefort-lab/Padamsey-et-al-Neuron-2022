import numpy as np

def check_setup_model_sim(modelparams, simparams):
	"""Run this before running loop simulations to check everything is in order"""

	from neuron import h
	from neuron.units import ms, mV
	import modelsetup as mset

	print('------------------------')
	print('Display model check...')

	# Create the  main compartment 
	# this works too:
	#soma, _, _ = modelsetup.setupsoma(modelparams, True)
	soma = mset.setupsoma(modelparams, True)[0]

	print("\t-Total duration of stimulation {:.2f} ms".format(simparams['stimduration']))
	print("\t-Total duration of simulation {:.2f} ms".format(simparams['simulation_length']))
	
	# Setup stim synapse:
	stimsyn = mset.Synapse(soma, simparams)

	# Check if inster syn double exp or if simulate it with play and Iclamp
	if simparams['intype']=='isyn':

		print("Use Iclamp for synapse input current")
		print(len(stimsyn.g_i_norm))

	#else:
		
	#print("\t-Reversal potential synaptic input {:.1f} mV".format(stimsyn.synE.e))

	stimsyn.set_ampmax(1) # dummy test

	# Check if artificial depo/hyperpo
	if modelparams.iclamprest:
		iclamprest = h.IClamp(soma(0.5))
		iclamprest.delay = 0 # ms
		iclamprest.dur = simparams['simulation_length']
		iclamprest.amp = modelparams.compute_currentclamprest(soma(0.5).area()) # uA/cm2 -> nA/cm2 -> nA # SA in um2 convert to cm2 # nA
		print("\t-Clamp resting potential {:.3f} uA/cm2".format(modelparams.iclamprest))

	print("\t-Syn noise for trial variability {:.3f} (sd)".format(stimsyn.syn_noise_sd))

	

def run_ginput_loop(inputs):
	"""Setup and run model synaptic input. For a given input, run N times to return results.
	Import here neuron module to allow multiprocessing (checked Neuron website Forum)."""

	from neuron import h
	from neuron.units import ms, mV
	import modelsetup as mset
    
    
	np.random.seed() # setup random seed (otherwise would use the same seed for every processor)
    
	h.dt = 0.01

	# Setup model and simulation parameters
	
	modelparams = inputs[0]
	simparams = inputs[1]
	ge_input = inputs[2]
	numloops = inputs[3]

	# Setup clipping AP threshold
	vm_clip = -40 # would be good to have this as param

	# Load the ... (see tuto)
	h.load_file("stdrun.hoc") # for run control

	# Create the  main compartment 
	soma = mset.setupsoma(modelparams)[0]

	# Prep the stim --------

	# Setup the synaptic input - excitatory only
	# Setup stim synapse:
	stimsyn = mset.Synapse(soma, simparams)
	#synE,stimE,ncE = mset.insertsyn(soma, simparams['synparams'])

	# Set the conductance conversion 
	#g_convert = soma(0.5).area()*1e-8 # uS  # SA in um2 convert to cm2

	# Check if artificial depo/hyperpo
	if modelparams.iclamprest:
		iclamprest = h.IClamp(soma(0.5))
		iclamprest.delay = 0 # ms
		iclamprest.dur = simparams['simulation_length']
		iclamprest.amp = modelparams.compute_currentclamprest(soma(0.5).area()) # uA/cm2 -> nA/cm2 -> nA # SA in um2 convert to cm2 # nA

	# Setup recording variables and AP threshold --------
	v_sim = h.Vector().record(soma(0.5)._ref_v) # Membrane potential vector
	t_sim = h.Vector().record(h._ref_t)  # Time stamp vector
	# Set the counter of action potentials
	apcounter = h.APCount(soma(0.5))
	# Change the threshold according to the model
	apcounter.thresh = modelparams.APthreshold
	# Set recording variables
	aptime_sim = h.Vector()  # AP times! count
	apcounter.record(aptime_sim)
  
	rec_var_run = {'apc':[], 'vmax':[], 'vh':[]}
	# Loop N times
	for iloop in range(numloops):

		# Set the syn conductances
		stimsyn.set_ampmax(ge_input)
		#ncE.weight[0] = ge_inputlist_convrt[iloop] # NetCon weight is a vector.

		# Init membrane potential
		h.finitialize(modelparams.E0 * mV)
		# Run
		h.continuerun(simparams['simulation_length'] * ms)

		# save apc
		rec_var_run['apc'].append(int(apcounter.n))

		# transform data to lists
		vm_list = list(v_sim) 
		t_list = list(t_sim)

		# save the max
		rec_var_run['vmax'].append(max(vm_list))
        
		# First remove everything above given threshold      
		vm_list_clipped = [vm for vm in vm_list if vm<vm_clip] 
		rec_var_run['vh'].append(np.quantile(vm_list_clipped, 0.975)) 

	return rec_var_run


def run_io_ginput(modelparams, simparams, list_input, numloops, doparallel=False, chunksize=1, max_workers=None):

	
	# Setup savelists
	rec_var =  {'apc':[], 'vmax':[], 'vh':[]}
	

	if doparallel:
		import concurrent.futures

		list_input_loop = [(modelparams, simparams, ge, numloops) for ge in list_input]

		with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
			results = list(executor.map(run_ginput_loop, list_input_loop, chunksize=chunksize))
		
			for key0 in rec_var.keys():
				rec_var[key0]=[result[key0] for result in results]

	else:
	
		for ge in list_input:
			rec_var_run = run_ginput_loop((modelparams, simparams, ge, numloops))
			
			# save all the new data points in the main list     
			for key0 in rec_var.keys():
				rec_var[key0].append(rec_var_run[key0])


	return rec_var


def get_firstspike_ginput(modelparams,simparams,offInput,onInput,precision):
	"""Finds the first g input, at a given precision, that generates a spike in given time window of simulation"""
    

	from neuron import h
	from neuron.units import ms, mV
	import modelsetup as mset

	h.dt = 0.01
    
	# Setup model and simulation parameters
	
	# Load the ... (see tuto)
	h.load_file("stdrun.hoc") # for run control

	# Create the  main compartment 
	soma = mset.setupsoma(modelparams)[0]

	# Prep the stim --------

	# Setup the synaptic input - excitatory only
	# Setup stim synapse:
	stimsyn = mset.Synapse(soma, simparams)
	#synE,stimE,ncE = mset.insertsyn(soma, simparams['synparams'])

	# Set the conductance conversion 
	#g_convert = soma(0.5).area()*1e-8 # uS  # SA in um2 convert to cm2

	# Check if artificial depo/hyperpo
	if modelparams.iclamprest:
		iclamprest = h.IClamp(soma(0.5))
		iclamprest.delay = 0 # ms
		iclamprest.dur = simparams['simulation_length']
		iclamprest.amp = modelparams.compute_currentclamprest(soma(0.5).area()) # uA/cm2 -> nA/cm2 -> nA # SA in um2 convert to cm2 # nA

	# Setup recording variables and AP threshold --------

	# Set the counter of action potentials
	apcounter = h.APCount(soma(0.5))
	# Change the threshold according to the model
	apcounter.thresh = modelparams.APthreshold
	# Set recording variables
	aptime_sim = h.Vector()  # AP times! count
	apcounter.record(aptime_sim)


	# Start the search for the first spike within given time window for given precision.
	#  We search by slicing the range of current in 2 until we find it.

	current_input = np.round(((offInput+onInput)/2)/(precision/10))*precision/10
	while (onInput-offInput)>precision :
		# Set the new test syn conductance
		stimsyn.set_ampmax(current_input)
		#ncE.weight[0] = current_input*g_convert # NetCon weight is a vector.

		# Init membrane potential
		h.finitialize(modelparams.E0 * mV)
		# Run
		h.continuerun(simparams['simulation_length'] * ms)

		# Check if has generated an AP
		if int(apcounter.n)==0:
			# If not AP, we need to increase the current
			offInput = current_input
			current_input = np.round(((current_input+onInput)/2)/(precision/10))*precision/10
			#current_input = (current_input+onInput)/2
		else:
			# If AP, we need to decrease the current
			onInput = current_input
			current_input = np.round(((current_input+offInput)/2)/(precision/10))*precision/10
			#current_input = (current_input+offInput)/2

	# Check if the on is on...
	stimsyn.set_ampmax(onInput)
	#ncE.weight[0] = onInput*g_convert # NetCon weight is a vector.
	# Run
	h.finitialize(modelparams.E0 * mV)
	h.continuerun(simparams['simulation_length'] * ms)
	n1 = int(apcounter.n)

	# ... and off is off 
	#ncE.weight[0] = offInput*g_convert # NetCon weight is a vector.
	stimsyn.set_ampmax(offInput)
	# Run
	h.finitialize(modelparams.E0 * mV)
	h.continuerun(simparams['simulation_length'] * ms)
	n0 = int(apcounter.n)

	is_success = n1>0 and n0==0

	# In the end, we return the 'on' input, and if success
	return onInput, offInput, is_success






def run_ginput_trace(modelparams, simparams, ge_input, output_g=False):

	from neuron import h
	from neuron.units import ms, mV
	import modelsetup as mset
    
	h.dt = 0.01

	# Setup model and simulation parameters

	# Load the ... (see tuto)
	h.load_file("stdrun.hoc") # for run control

	# Create the  main compartment 
	soma = mset.setupsoma(modelparams)[0]

	# Prep the stim --------

	# Setup the synaptic input - excitatory only
	# Setup stim synapse:
	stimsyn = mset.Synapse(soma, simparams)
	#synE,stimE,ncE = mset.insertsyn(soma, simparams['synparams'])

	# Set the conductance conversion 
	#g_convert = soma(0.5).area()*1e-8 # uS  # SA in um2 convert to cm2

	# Check if artificial depo/hyperpo
	if modelparams.iclamprest:
		iclamprest = h.IClamp(soma(0.5))
		iclamprest.delay = 0 # ms
		iclamprest.dur = simparams['simulation_length']
		iclamprest.amp = modelparams.compute_currentclamprest(soma(0.5).area()) # uA/cm2 -> nA/cm2 -> nA # SA in um2 convert to cm2 # nA

	# Setup recording variables and AP threshold --------
	v_sim = h.Vector().record(soma(0.5)._ref_v) # Membrane potential vector
	t_sim = h.Vector().record(h._ref_t)  # Time stamp vector
	# Set the counter of action potentials
	apcounter = h.APCount(soma(0.5))
	# Change the threshold according to the model
	apcounter.thresh = modelparams.APthreshold
	# Set recording variables
	aptime_sim = h.Vector()  # AP times! count
	apcounter.record(aptime_sim)
    
	if output_g:
		synE_g_sim = h.Vector().record(stimsyn.synE._ref_g)
  
	# Prep the syn conductances
	stimsyn.set_ampmax(ge_input)
	#ncE.weight[0] = ge_input*g_convert # NetCon weight is a vector.

	# Setup savelists
	if output_g:
		trace_var =  {'t':[], 'apc':[], 'v':[], 'syn_ge':[]} 
	else:
		trace_var =  {'t':[], 'apc':[], 'v':[]} 

	
	# Init membrane potential
	h.finitialize(modelparams.E0 * mV)
	# Run
	h.continuerun(simparams['simulation_length'] * ms)

	t_list = list(t_sim)
	vm_list = list(v_sim)
    
	# save apc
	trace_var['apc'].append(int(apcounter.n))
	trace_var['v'].append(vm_list)
    
	if output_g:
		trace_var['syn_ge'].append(list(synE_g_sim))

	trace_var['t'] = list(t_sim)

	return trace_var