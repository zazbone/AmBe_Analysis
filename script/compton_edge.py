import numpy as np
import matplotlib.pyplot as plt

me = 0.511  # MeV
EI = 4.438  # MeV

def calc_EF(EI, theta):
    num = EI * me
    den = 1 - np.cos(theta)
    den *= EI
    den += me
    return num / den


theta = np.linspace(0, np.pi, 255)
EF = calc_EF(EI, theta)
print(np.min(EF), np.max(EF))
plt.plot(theta, EF, color="red", label=f"$\\gamma$ energy (Min: {np.min(EF):.2E})")
plt.plot(theta, EI- EF, color="blue", label=f"$e^-$ energy (Min: {np.max(EI - EF):.2E})")
plt.hlines(0, 0, np.pi, color="black", linestyle="dashed")
plt.hlines(EI, 0, np.pi, color="black", linestyle="dashed")
plt.text(-0.41, EI, "$E_{\gamma I}$=4.438 --")
plt.xlabel(r"$\theta$ scattering angle (rd)")
plt.ylabel(r" Final energy (MeV)")
plt.legend(loc="best")
plt.show()