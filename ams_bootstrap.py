"""
Moyenne tensorielle par bootstrap : le Fortran d'origine (`bootams`,
AMS_OSX_x.f95:391-393) est un STUB VIDE - `subroutine bootams; return;
end` ne fait litteralement rien ("programme en construction" comme
`amsgmt`/`grandcerc`). Il n'y a donc rien a porter depuis le Fortran ici.

A la demande explicite de l'utilisateur, recherche faite dans la
distribution PyPI de PmagPy (`pip install pmagpy`, verifie disponible et
installe - voir `pmagpy.pmag`) : `s_boot`/`sbootpars` implementent
exactement ce que le menu "ellipses Bootstrap" aurait du faire - bootstrap
non parametrique (par defaut) ou parametrique sur une liste de tenseurs
6-composantes (meme convention k11,k22,k33,k12,k23,k13 que le format .ANI
AMS), avec ajustement d'une distribution de Kent sur les vecteurs propres
rechantillonnes pour obtenir une ellipse de confiance (zeta/eta = demi-
grand/demi-petit axe en degres, chacun avec sa propre orientation
dec/inc) - une methode DIFFERENTE de l'approche Jelinek (1978) deja portee
(`ams_stats.ellips`), mais un standard reconnu en paleomagnetisme
(Constable & Tauxe 1990) pour le meme probleme.

Ce module reutilise directement `pmagpy.pmag` (paquet installe, pas une
reimplementation) - contrairement au reste du port AMS_Py qui retranscrit
le Fortran ligne a ligne, il n'y a ici aucune fidelite a preserver
puisque la routine Fortran d'origine n'a jamais rien fait."""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from pmagpy import pmag

from ams_selection import AMSMeasurement, apply_orientation


@dataclass
class BootstrapAxisResult:
    dec: float
    inc: float
    eigenvalue_mean: float
    eigenvalue_sigma: float
    zeta: float       # demi-grand axe de l'ellipse de confiance (deg)
    zeta_dec: float
    zeta_inc: float
    eta: float         # demi-petit axe (deg)
    eta_dec: float
    eta_inc: float
    draws: List[Tuple[float, float]] = field(default_factory=list)  # (dec,inc) de chaque tirage bootstrap - equivalent de 'ivec=True' cote PmagPy (plot_aniso)


@dataclass
class BootstrapResult:
    n: int = 0
    nb: int = 0
    parametric: bool = False
    axes: List[BootstrapAxisResult] = field(default_factory=list)  # [k1, k2, k3]
    id: str = ""


def compute_bootstrap_mean(
    measurements: List[AMSMeasurement],
    orientation: int,
    nb: int = 1000,
    parametric: bool = False,
    random_seed=None,
) -> Optional[BootstrapResult]:
    """Equivalent de `pmag.s_boot`+`pmag.sbootpars` applique a la liste de
    mesures `measurements` (orientees selon `orientation`, meme convention
    que `ams_stats.tsmean`). `nb` : nombre de tirages bootstrap (defaut
    1000, comme le Fortran `bootstrap` du menu Calcul de Starmac). Renvoie
    None si moins de 3 mesures (meme garde que `tsmean`)."""
    if len(measurements) < 3:
        return None

    Ss = []
    for m in measurements:
        a = apply_orientation(m, orientation)
        Ss.append([a[0, 0], a[1, 1], a[2, 2], a[0, 1], a[1, 2], a[0, 2]])

    Tmean, Vmean, Taus, Vs = pmag.s_boot(Ss, ipar=1 if parametric else 0, nb=nb, random_seed=random_seed)
    bpars = pmag.sbootpars(Taus, Vs)

    axes = []
    for i in range(3):
        tau = Tmean[i]
        tau = float(tau.real) if hasattr(tau, "real") else float(tau)
        draws = [(float(draw[i][0]), float(draw[i][1])) for draw in Vs]
        axes.append(BootstrapAxisResult(
            dec=float(Vmean[i][0]), inc=float(Vmean[i][1]),
            eigenvalue_mean=tau,
            eigenvalue_sigma=float(bpars[f"t{i + 1}_sigma"]),
            zeta=float(bpars[f"v{i + 1}_zeta"]), zeta_dec=float(bpars[f"v{i + 1}_zeta_dec"]),
            zeta_inc=float(bpars[f"v{i + 1}_zeta_inc"]),
            eta=float(bpars[f"v{i + 1}_eta"]), eta_dec=float(bpars[f"v{i + 1}_eta_dec"]),
            eta_inc=float(bpars[f"v{i + 1}_eta_inc"]),
            draws=draws,
        ))
    return BootstrapResult(n=len(measurements), nb=nb, parametric=parametric, axes=axes)


def format_bootstrap_result(res: BootstrapResult) -> str:
    method = "parametrique" if res.parametric else "non parametrique"
    lines = [
        f" Bootstrap tensorial mean (PmagPy pmag.s_boot/sbootpars, {method}, "
        f"nb={res.nb}, n={res.n})",
    ]
    for i, ax in enumerate(res.axes, start=1):
        lines.append(
            f" k{i}: tau={ax.eigenvalue_mean:.5f} (sigma={ax.eigenvalue_sigma:.5f})   "
            f"dec={ax.dec:6.1f}  inc={ax.inc:5.1f}   "
            f"zeta={ax.zeta:5.1f} (dec={ax.zeta_dec:6.1f} inc={ax.zeta_inc:5.1f})   "
            f"eta={ax.eta:5.1f} (dec={ax.eta_dec:6.1f} inc={ax.eta_inc:5.1f})"
        )
    return "\n".join(lines) + "\n"
