import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Parameters
# -------------------------
frh_crit = 70.
cap_max_in = 150.0
cap_max_inc = 90.0

frh = np.linspace(0.0, 1.2, 500)

# -------------------------
# fx formulations
# -------------------------
def fx_linear(frh):
    fx = 2.0 * (frh - frh_crit)
    return np.clip(fx, -1.0, 1.0)

def fx_quadratic(frh):
    fx = 4.0 * (frh - frh_crit) * np.abs(frh - frh_crit)
    return np.clip(fx, -1.0, 1.0)

def fx_exponential(frh):
    fx = (2.0 / 0.78) * np.exp(-(frh - frh_crit)**2) * (frh - frh_crit)
    return np.clip(fx, -1.0, 1.0)

def fx_tanh(frh):
    fx = np.tanh(2.0 * (frh - frh_crit))
    return np.clip(fx, -1.0, 1.0)

# -------------------------
# cap_max computation
# -------------------------
def compute_cap_max(fx):
    cap = cap_max_in + fx * cap_max_inc
    return np.clip(cap, 10.0, 150.0)

# Compute fx
fx_lin = fx_linear(frh)
fx_quad = fx_quadratic(frh)
fx_exp = fx_exponential(frh)
fx_t = fx_tanh(frh)

# Compute cap_max
cap_lin = compute_cap_max(fx_lin)
cap_quad = compute_cap_max(fx_quad)
cap_exp = compute_cap_max(fx_exp)
cap_t = compute_cap_max(fx_t)

# -------------------------
# Plot
# -------------------------
plt.figure(figsize=(8, 6))

plt.plot(frh, cap_lin, label='Linear')
plt.plot(frh, cap_quad, label='Quadratic')
plt.plot(frh, cap_exp, label='Exponential')
plt.plot(frh, cap_t, label='Tanh')

# Reference lines
plt.axvline(frh_crit, linestyle='--', label='frh_crit')
plt.axhline(10, linestyle=':', linewidth=1)
plt.axhline(150, linestyle=':', linewidth=1)

# Labels
plt.xlabel('frh')
plt.ylabel('cap_max')
plt.title('cap_max vs frh for Different fx Formulations')
plt.legend()
plt.grid()

plt.show()
