"""Mass database of nuclei, single-Lambda and double-Lambda hypernuclei.

All masses are nuclear masses in MeV/c^2.
"""

from dataclasses import dataclass

# Physical constants (MeV/c^2), PDG / CODATA
AMU = 931.49410242
M_ELECTRON = 0.51099895
M_NEUTRON = 939.56542052
M_PROTON = 938.27208816
M_LAMBDA = 1115.683
M_XI_MINUS = 1321.71
M_PI_MINUS = 139.57039
M_PI_ZERO = 134.9768

# Bond energy of two Lambdas in the s-shell (Nagara event, MeV)
DELTA_B_LL = 0.67

ELEMENTS = ["n", "H", "He", "Li", "Be", "B", "C", "N", "O"]

# Atomic mass excess (keV), AME2020. (A, Z): (excess, particle-bound)
# Unbound nuclei are kept only as cores of hypernuclei.
MASS_EXCESS_KEV = {
  (1, 0): (8071.318, True),
  (1, 1): (7288.971, True),
  (2, 1): (13135.723, True),
  (3, 1): (14949.811, True),
  (3, 2): (14931.219, True),
  (4, 2): (2424.916, True),
  (5, 2): (11231.0, False),
  (5, 3): (11679.0, False),
  (6, 2): (17592.10, True),
  (6, 3): (14086.88, True),
  (6, 4): (18375.0, False),
  (7, 2): (26073.0, False),
  (7, 3): (14907.10, True),
  (7, 4): (15769.00, True),
  (8, 2): (31609.7, True),
  (8, 3): (20945.80, True),
  (8, 4): (4941.67, False),
  (8, 5): (22921.6, True),
  (9, 3): (24954.9, True),
  (9, 4): (11348.45, True),
  (9, 5): (12416.5, False),
  (9, 6): (28910.5, True),
  (10, 4): (12607.49, True),
  (10, 5): (12050.61, True),
  (10, 6): (15698.6, True),
  (11, 3): (40728.0, True),
  (11, 4): (20177.2, True),
  (11, 5): (8667.71, True),
  (11, 6): (10650.3, True),
  (12, 4): (25077.8, True),
  (12, 5): (13369.4, True),
  (12, 6): (0.0, True),
  (12, 7): (17338.1, True),
  (13, 5): (16562.1, True),
  (13, 6): (3125.009, True),
  (13, 7): (5345.48, True),
  (14, 4): (39950.0, True),
  (14, 5): (23664.0, True),
  (14, 6): (3019.89, True),
  (14, 7): (2863.417, True),
  (15, 5): (28958.0, True),
  (15, 6): (9873.1, True),
  (15, 7): (101.439, True),
  (16, 6): (13694.0, True),
  (16, 7): (5683.9, True),
  (16, 8): (-4737.00, True),
}

# Measured B_Lambda (MeV) of ground states, keyed by hypernucleus (A, Z).
# Emulsion compilation (Juric 1973 / Davis 2005) unless noted.
B_LAMBDA_MEASURED = {
  (3, 1): 0.13,
  (4, 1): 2.04,
  (4, 2): 2.39,
  (5, 2): 3.12,
  (6, 2): 4.18,
  (7, 2): 5.68,   # JLab (e,e'K+)
  (8, 2): 7.16,
  (7, 3): 5.58,
  (8, 3): 6.80,
  (9, 3): 8.50,
  (7, 4): 5.16,
  (8, 4): 6.84,
  (9, 4): 6.71,
  (10, 4): 9.11,
  (9, 5): 8.29,
  (10, 5): 8.89,
  (11, 5): 10.24,
  (12, 5): 11.37,
  (12, 6): 10.76,
  (13, 6): 11.69,
  (14, 6): 12.17,
  (15, 7): 13.59,
  (16, 7): 13.76,  # JLab (e,e'K+)
}

# Measured B_LambdaLambda (MeV), keyed by (A, Z).
B_LL_MEASURED = {
  (6, 2): 6.91,  # Nagara
}

MAX_Z = 7
MAX_A = 17


@dataclass(frozen=True)
class Species:
  name: str
  a: int        # baryon number
  z: int        # charge
  n_lambda: int
  mass: float   # MeV/c^2
  estimated: bool = False


def nuclear_mass(a, z):
  excess = MASS_EXCESS_KEV[(a, z)][0] * 1e-3
  return a * AMU + excess - z * M_ELECTRON


def nucleus_name(a, z):
  if (a, z) == (1, 0):
    return "n"
  if (a, z) == (1, 1):
    return "p"
  if (a, z) == (2, 1):
    return "d"
  if (a, z) == (3, 1):
    return "t"
  return f"{a}{ELEMENTS[z]}"


def hyper_name(a, z, n_lambda):
  tag = "L" * n_lambda
  return f"{tag}{a}{ELEMENTS[z]}"


def isobar_mean_b_lambda():
  """Mean measured B_Lambda for each A, linearly interpolated."""
  by_a = {}
  for (a, _), b in B_LAMBDA_MEASURED.items():
    by_a.setdefault(a, []).append(b)
  mean = {a: sum(v) / len(v) for a, v in by_a.items()}
  known = sorted(mean)
  for a in range(known[0], known[-1] + 1):
    if a in mean:
      continue
    lo = max(k for k in known if k < a)
    hi = min(k for k in known if k > a)
    frac = (a - lo) / (hi - lo)
    mean[a] = mean[lo] + frac * (mean[hi] - mean[lo])
  return mean


def build_candidates(use_estimates=True):
  """Return nuclei, single-L and double-L candidates (not yet checked)."""
  nuclei = [
    Species(nucleus_name(a, z), a, z, 0, nuclear_mass(a, z))
    for (a, z), (_, bound) in MASS_EXCESS_KEV.items()
    if bound and z <= MAX_Z
  ]
  lam = Species("L", 1, 0, 1, M_LAMBDA)

  b_est = isobar_mean_b_lambda()
  singles = {}
  for (ac, z) in MASS_EXCESS_KEV:
    a = ac + 1
    if z < 1 or z > MAX_Z or a > MAX_A:
      continue
    measured = (a, z) in B_LAMBDA_MEASURED
    if not measured and not (use_estimates and a in b_est):
      continue
    b = B_LAMBDA_MEASURED[(a, z)] if measured else b_est[a]
    mass = nuclear_mass(ac, z) + M_LAMBDA - b
    singles[(a, z)] = (Species(hyper_name(a, z, 1), a, z, 1, mass,
                               not measured), b)

  doubles = []
  for (ac, z) in MASS_EXCESS_KEV:
    a = ac + 2
    if z < 1 or z > MAX_Z or a > MAX_A:
      continue
    if (a, z) in B_LL_MEASURED:
      b_ll, est = B_LL_MEASURED[(a, z)], False
    elif (a - 1, z) in singles:
      sp, b_l = singles[(a - 1, z)]
      b_ll, est = 2.0 * b_l + DELTA_B_LL, True
    else:
      continue
    mass = nuclear_mass(ac, z) + 2.0 * M_LAMBDA - b_ll
    doubles.append(Species(hyper_name(a, z, 2), a, z, 2, mass, est))

  return nuclei, lam, [s for s, _ in singles.values()], doubles
