# Noise-induced-firing-patterns-in-HH-neuron⚡

Lets first understand what does noise induced firing patterns and Hodgkin-Huxley neuron model mean.<br>
This project studies how noise influences the firing pattern of a Hodgkin-Huxley (HH) neuron 
- **Noise-induced firing patterns:**<br>
   Random fluctuations in the neuron's input can trigger or modify the timing and frequency of the action potentioal.
- **Hodgkin-Huxley (HH) model:**<br>
   A mathematical model describing how neurons generate action potentials using voltage-dependet ion channels(Na+, K+, leak currents)

# Project Extensions

The project extends the classical HH model by including:<br>
- **Stochastic input (Noise)**<br>
    Gaussian noise is added to the input current to study its effects on spike         timing.<br>
- **Calcium (Ca²⁺) cuurents**<br>
   Voltage-gated calcium channels that contribute to depolarization.<br>
- **Calcium-activated potssium (KCa) currents**<br>
   Potassium channels that open in response to intracellular calcium, aiding          repolarization.<br>
- **Intracellular calcium dynamics:**<br>
   Simulates changes in intraceelular calcium concentration over time.<br>
- **Multiple trials:**<br>
   The simulation can run multiple trials with noise to visualize variability in      firing patterns.<br>

## Model Outputs

The model generates:

- Action potential plots ( with and without Ca²⁺/KCacurrents)
- Gating variable plots  (m, h, n, p)
- Ionic cuurent plots (INa, IK, ICa, IKCa, IL)
- Intracellular calcium dynamics
- Noise-induced spike train plots
- Firing rate in Hz
- 3D plot

# Files and Contents

- **Firing_patterns.py**<br>
    It contains the python code implementing the HH neuron model with noise, calcium and KCa currents. It also includes functionality of multiple trials,       spike detection , ISI calculation and plots.
  
- **Plots**<br>
    It contains ISI histograms AND 3D plot showing the relationship between time, membrane potential and calcium activation
