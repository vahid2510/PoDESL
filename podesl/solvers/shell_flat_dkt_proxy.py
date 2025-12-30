import numpy as np
from .plate_mindlin_winkler import solve_plate_mindlin_winkler as _pm
# Thin flat shell proxy: route to plate mindlin without foundation (kw=0 by default)
def solve_shell_flat_dkt_static(env):
    env2 = dict(env)
    env2.setdefault("kw", 0.0)  # no foundation
    return _pm(env2)
