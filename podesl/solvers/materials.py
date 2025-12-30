from __future__ import annotations
import numpy as np
from dataclasses import dataclass

@dataclass
class VonMises1D:
    E: float
    nu: float
    sigma_y0: float
    H: float
    # state
    eps_p: float = 0.0
    alpha: float = 0.0   # equivalent plastic strain
    sigma: float = 0.0

    def step(self, deps: float):
        E = self.E; H = self.H; sigy0 = self.sigma_y0
        sigma_trial = self.sigma + E*deps
        f = abs(sigma_trial) - (sigy0 + H*self.alpha)
        if f <= 0.0:
            self.sigma = sigma_trial
            Et = E
            return self.sigma, Et, {"elastic": True, "eps_p": self.eps_p, "alpha": self.alpha}
        dgamma = f / (E + H)
        sign = 1.0 if sigma_trial >= 0.0 else -1.0
        self.sigma = sigma_trial - dgamma*E*sign
        self.eps_p += dgamma*sign
        self.alpha += dgamma
        Et = (E*H)/(E+H)
        return self.sigma, Et, {"elastic": False, "eps_p": self.eps_p, "alpha": self.alpha}

def uniaxial_von_mises_curve(env):
    E = float(env["E"]); nu = float(env.get("nu", 0.3))
    sigma_y0 = float(env.get("sigma_y0", env.get("sigy0", 250e6)))
    H = float(env.get("H", 1.0e9))
    eps_max = float(env.get("eps_max", 0.02))
    nsteps = int(env.get("nsteps", 200))
    path = env.get("eps_path", None)
    if path is None:
        eps = np.linspace(0.0, eps_max, nsteps)
    else:
        eps = np.array(path, dtype=float).reshape(-1)
        nsteps = eps.size
    mat = VonMises1D(E=E, nu=nu, sigma_y0=sigma_y0, H=H)
    sigma = np.zeros(nsteps)
    eps_p = np.zeros(nsteps)
    Et = np.zeros(nsteps)
    prev = 0.0
    for i,e in enumerate(eps):
        deps = e - prev
        s, et, st = mat.step(deps)
        sigma[i] = s; eps_p[i] = st["eps_p"]; Et[i] = et
        prev = e
    return {"eps": eps, "sigma": sigma, "eps_p": eps_p, "Et": Et}
