import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs("plots", exist_ok=True)

dt = 0.01          # time step (ms)
t_max = 50         # total time (ms)
t = np.arange(0, t_max, dt)

# Membrane capacitance (µF/cm^2)
Cm = 1.0

# Maximum conductances (mS/cm^2)
gNa = 120.0   # Sodium
gK = 36.0     # Potassium
gL = 0.3      # Leak
gCa = 1.0     #calcium
gKCa = 5.0  # Calcium-activated K+ current


# Reversal potentials (mV)
ENa = 50.0
EK = -77.0
EL = -54.4
ECa = 120.0

# Calcium dynamics parameters
alpha_Ca = 0.0001  # converts I_Ca to concentration
tau_Ca = 80.0      # ms, decay of intracellular Ca
Kd_Ca = 0.5        # µM, half-activation for KCa
n_Hill = 2         # Hill coefficient

I_ext = np.zeros(len(t))
I_ext[(t > 10) & (t < 40)] = 10  # µA/cm² stimulus

sigma = 2.0  # noise amplitude
I_noisy = I_ext + sigma * np.random.randn(len(I_ext))  # add Gaussian noise

## sodium

def alpha_m(V):
    return (0.1 * (V + 40)) / (1 - np.exp(-(V + 40) / 10))

def beta_m(V):
    return 4 * np.exp(-(V + 65) / 18)

## potassium

def alpha_h(V):
    return 0.07 * np.exp(-(V + 65) / 20)

def beta_h(V):
    return 1 / (1 + np.exp(-(V + 35) / 10))

## leak

def alpha_n(V):
    return (0.01 * (V + 55)) / (1 - np.exp(-(V + 55) / 10))

def beta_n(V):
    return 0.125 * np.exp(-(V + 65) / 80)

## calcium

def alpha_p(V):
    return 0.1 * (V + 25) / (1 - np.exp(-(V + 25) / 10))

def beta_p(V):
    return 0.1 * np.exp(-(V + 35) / 10)

def run_HH_Ca_KCa(include_Ca=True, include_KCa= True, I_input=None):
    if I_input is None:
        I_input = I_ext
    V = -65 * np.ones(len(t))

    m = np.zeros(len(t))
    h = np.zeros(len(t))
    n = np.zeros(len(t))
    p = np.zeros(len(t))
    Ca_i = np.zeros(len(t))

    m[0] = alpha_m(V[0]) / (alpha_m(V[0]) + beta_m(V[0]))
    h[0] = alpha_h(V[0]) / (alpha_h(V[0]) + beta_h(V[0]))
    n[0] = alpha_n(V[0]) / (alpha_n(V[0]) + beta_n(V[0]))
    p[0] = alpha_p(V[0]) / (alpha_p(V[0]) + beta_p(V[0]))
    Ca_i[0] = 0.0  # initial intracellular Ca

    INa_trace = np.zeros(len(t))
    IK_trace  = np.zeros(len(t))
    IL_trace  = np.zeros(len(t))
    ICa_trace = np.zeros(len(t))
    IKCa_trace = np.zeros(len(t))

    for i in range(1, len(t)):
        m[i] = m[i-1] + dt * (alpha_m(V[i-1])*(1-m[i-1]) - beta_m(V[i-1])*m[i-1])
        h[i] = h[i-1] + dt * (alpha_h(V[i-1])*(1-h[i-1]) - beta_h(V[i-1])*h[i-1])
        n[i] = n[i-1] + dt * (alpha_n(V[i-1])*(1-n[i-1]) - beta_n(V[i-1])*n[i-1])
        p[i] = p[i-1] + dt * (alpha_p(V[i-1])*(1-p[i-1]) - beta_p(V[i-1])*p[i-1])

        INa = gNa * (m[i-1]**3) * h[i-1] * (V[i-1] - ENa)
        IK  = gK  * (n[i-1]**4) * (V[i-1] - EK)
        IL  = gL  * (V[i-1] - EL)
        ICa = gCa * (p[i-1]**2) * (V[i-1] - ECa) if include_Ca else 0
        IKCa = gKCa*(Ca_i[i-1]**n_Hill/(Ca_i[i-1]**n_Hill + Kd_Ca**n_Hill))*(V[i-1]-EK) if include_KCa else 0


       # Update intracellular calcium
        dCa = -alpha_Ca*ICa - Ca_i[i-1]/tau_Ca
        Ca_i[i] = Ca_i[i-1] + dt*dCa

        INa_trace[i] = INa
        IK_trace[i]  = IK
        IL_trace[i]  = IL
        ICa_trace[i] = ICa
        IKCa_trace[i] = IKCa

        V[i] = V[i-1] + dt * (I_input[i] - INa - IK - IL - ICa- IKCa) / Cm

    return V, m, h, n, p, INa_trace, IK_trace, IL_trace, ICa_trace, IKCa_trace, Ca_i

# --- Multiple trials with noise ---
num_trials = 5
plt.figure(figsize=(10,4))
for trial in range(num_trials):
    V_trial, *_ = run_HH_Ca_KCa(include_Ca=True, include_KCa=True, I_input=I_noisy)
    plt.plot(t, V_trial, label=f'Trial {trial+1}')
plt.xlabel("Time (ms)")
plt.ylabel("Membrane Potential (mV)")
plt.title("HH + Ca/KCa with Noise")
plt.legend()
plt.tight_layout()
plt.show()

V_to_analyze = V_trial
threshold = 0

spike_times = t[1:][(V_to_analyze[:-1] < threshold) & (V_to_analyze[1:] >= threshold)]

ISI = np.diff(spike_times)

plt.figure()
plt.hist(ISI, bins=10)
plt.xlabel("Inter-Spike Interval (ms)")
plt.ylabel("Count")
plt.title("ISI Histogram")
plt.show()

firing_rate = len(spike_times) / (t_max/1000)
print("Firing rate (Hz):", firing_rate)

# Run model WITHOUT calcium
V_noCa, *_ = run_HH_Ca_KCa(include_Ca=False, include_KCa=False)

# Run model WITH calcium
V, m, h, n, p, INa_trace, IK_trace, IL_trace, ICa_trace, IKCa_trace, Ca_i = run_HH_Ca_KCa(include_Ca=True, include_KCa=True)

plt.figure(figsize=(10,4))
plt.plot(t, V_noCa, label="Without Ca²⁺")
plt.plot(t, V, label="With Ca²⁺ + KCa")
plt.xlabel("Time (ms)")
plt.ylabel("Membrane Potential (mV)")
plt.title("Comparison of Action Potentials")
plt.legend()
plt.tight_layout()
plt.savefig("plots/action_potential_comparison.png", dpi=300)
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(t, V)
plt.xlabel("Time (ms)")
plt.ylabel("Membrane Potential (mV)")
plt.title("HH + Ca²⁺ + KCa Action Potential")
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(t, m, label="m (Na activation)")
plt.plot(t, h, label="h (Na inactivation)")
plt.plot(t, n, label="n (K activation)")
plt.plot(t, p, label="p (Ca activation)")
plt.xlabel("Time (ms)")
plt.ylabel("Gating Value")
plt.legend()
plt.title("Gating Variables")
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(t, INa_trace, label="INa")
plt.plot(t, IK_trace, label="IK")
plt.plot(t, ICa_trace, label="ICa")
plt.plot(t, IKCa_trace, label="IKCa")
plt.plot(t, IL_trace, label="IL")
plt.xlabel("Time (ms)")
plt.ylabel("Current (µA/cm²)")
plt.legend()
plt.title("Ionic Currents")
plt.tight_layout()
plt.savefig("plots/ionic_currents.png", dpi=300)
plt.show()

plt.figure(figsize=(10,4))
plt.plot(t, Ca_i)
plt.xlabel("Time (ms)")
plt.ylabel("Intracellular Ca²⁺ (µM)")
plt.title("Intracellular Calcium Dynamics")
plt.tight_layout()
plt.savefig("plots/intracellular_calcium.png", dpi=300)
plt.show()

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

ax.plot(t, V, p, linewidth=2)

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Membrane Potential (mV)')
ax.set_zlabel('Ca Activation (p)')
ax.set_title('3D Relationship Between Time, Voltage and Calcium Activation')

plt.tight_layout()
plt.savefig("plots/3d_voltage_calcium.png", dpi=300)
plt.show()

