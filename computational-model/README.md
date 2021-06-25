# Computational model for Padamsey et al.

Neuron model: Hodgkin-Huxley type, single compartment.

Numerical simulations were performed using the NEURON simulation environment interfaced with Python (https://www.neuron.yale.edu/neuron/). Analysis was done in Python, channel calibration in MATLAB.


## Simulations

We ran the simulations and visualised results with jupyter notebooks (.ipynb files). The notebooks are in the folder [neuron-model](neuron-model). Note that clicking on the .ipynb files via GitHub gives a preview.

| File           | Description | 
| ---------------    | ----------- | 
| `simu_toy_model_passive.ipynb`  |  Toy model (Figure 3/ H-K) <br />Passive components only, no voltage-gated channels |
| `simu_deterministic_model.ipynb`  |  Deterministic channels model  (Figure 5/ C-F) <br />- Simulate synaptic input trial-trial variability <br />- Input/Output and tuning curves|
| `simu_stochastic_model.ipynb`  |  Stochastic channels model  (Figure 5/ I-L) <br />- Voltage-dependent membrane noise <br />- Input/Output and tuning curves|
| `simu_trace_detchannels.ipynb`  |  Deterministic channels model  (Figure 5/ B) <br />- Simulate synaptic input trial-trial variability <br />- Temporal traces|
| `simu_trace_stochannels.ipynb`  |  Stochastic channels model  (Figure 5/ I-L) <br />- Voltage-dependent membrane noise <br />- Temporal traces|


## Channels models

Custom channels were defined using NEURON's Channel Builder:

https://www.neuron.yale.edu/neuron/static/docs/chanlbild/main.html

All channels are saved in a single session file, and can be found in the folder [neuron-model](neuron-model).

| File           | Description | 
| ---------------    | ----------- | 
| `ch_carter_subtchan.ses`  |  Deterministic channels |
| `stch_carter_subtchan.ses`  |  Stochastic channels |


