# Computational model for Padamsey et al.

Neuron model: Hodgkin-Huxley type, single compartment.

Numerical simulations were performed using the NEURON simulation environment interfaced with Python (https://www.neuron.yale.edu/neuron/). Analysis was done in Python, channel calibration in MATLAB.


## Simulations

We ran the simulations and visualised results with jupyter notebooks (.ipynb files). The notebooks are in the folder [neuron-model](neuron-model). Note that clicking on the .ipynb files via GitHub gives a preview.

| File           | Description | 
| ---------------    | ----------- | 
| `simu_toy_model_passive.ipynb`  |  Toy model (Figure 3/ H-K) <br />Passive components only, no voltage-gated channels |
| `simu_deterministic_model.ipynb`  |  Input/Output and tuning curves (Figure 5/ C-F) <br />- Deterministic channels model <br />- Simulate synaptic input trial-trial variability|
| `simu_stochastic_model.ipynb`  |  Input/Output and tuning curves (Figure 5/ I-L) <br />- Stochastic channels model<br />- Voltage-dependent membrane noise|
| `simu_trace_detchannels.ipynb`  |  Temporal traces (Figure 5/ B) <br />- Deterministic channels model <br />- Simulate synaptic input trial-trial variability|
| `simu_trace_stochannels.ipynb`  |  Temporal traces (Figure 5/ H) <br />- Stochastic channels model <br />- Voltage-dependent membrane noise |


## Channels models

Custom channels were defined using NEURON's Channel Builder:

https://www.neuron.yale.edu/neuron/static/docs/chanlbild/main.html

All channels are saved in a single session file, and can be found in the folder [neuron-model](neuron-model).

| File           | Description | 
| ---------------    | ----------- | 
| `ch_carter_subtchan.ses`  |  Deterministic channels |
| `stch_carter_subtchan.ses`  |  Stochastic channels |


