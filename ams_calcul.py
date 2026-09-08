"""
Calculs AMS : port des routines de `Util_AMS.f95` qui ne dependent PAS des
statistiques tensorielles de `anisotropie.f` (tsmean/Hext, en cours de
portage separement). Pour l'instant : `corpaldir` (correction inverse d'une
direction paleomagnetique par un tenseur AMS) et `soustract` (soustraction de
deux tenseurs, protocole de type ARM/susceptibilite anhysteretique).
"""

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

from ams_selection import AMSMeasurement, cart, polere


def inverse_symmetric_3x3(
    a: float, b: float, c: float, d: float, e: float, f: float,
    x: float, y: float, z: float,
) -> Tuple[float, float, float]:
    """Equivalent de `inverse` (Util_AMS.f95:120-136) : applique l'inverse de
    la matrice symetrique [[a,b,c],[b,d,e],[c,e,f]] au vecteur (x,y,z).
    Formule identique (verifiee terme a terme) a `calcul.inverse_symmetric_3x3`
    du port STARpaleomag_Py - meme heritage de code source."""
    det = a * d * f - a * e * e - b * b * f + 2 * b * c * e - c * c * d
    if det == 0:
        return x, y, z
    a11, a12, a13 = (d * f - e * e) / det, -(b * f - c * e) / det, (b * e - c * d) / det
    a22, a23 = (a * f - c * c) / det, -(a * e - b * c) / det
    a33 = (a * d - b * b) / det
    x1 = a11 * x + a12 * y + a13 * z
    y1 = a12 * x + a22 * y + a23 * z
    z1 = a13 * x + a23 * y + a33 * z
    return x1, y1, z1


def correct_direction_with_tensor(
    tensor: AMSMeasurement, intensity: float, dec: float, inc: float,
) -> Tuple[float, float, float]:
    """Equivalent du coeur de `corpaldir` (Util_AMS.f95:2-115) : corrige une
    direction paleomagnetique (intensite,dec,inc) de l'anisotropie donnee par
    `tensor` (le PREMIER tenseur de la liste selectionnee, comme le Fortran -
    "attention prend le premier tenseur de la liste", pas une moyenne).
    Normalise le tenseur par sa trace/3 avant d'appliquer l'inverse (meme
    convention que `calcul.compute_anicor_factor` cote STARpaleomag_Py), puis
    reconvertit en (intensite, dec, inc) via `polere`."""
    x, y, z = cart(intensity, dec, inc)
    rnorm = (tensor.k11 + tensor.k22 + tensor.k33) / 3.0
    a1 = tensor.k11 / rnorm
    b1 = tensor.k12 / rnorm
    c1 = tensor.k13 / rnorm
    d1 = tensor.k22 / rnorm
    e1 = tensor.k23 / rnorm
    f1 = tensor.k33 / rnorm
    x1, y1, z1 = inverse_symmetric_3x3(a1, b1, c1, d1, e1, f1, x, y, z)
    return polere(x1, y1, z1)


@dataclass
class SubtractionResult:
    """Les 3 tenseurs derives de `soustract` (Util_AMS.f95:158-283) pour
    exactement 2 mesures selectionnees - PAS la variante liste
    (`soustractliste`, plusieurs paires consecutives de meme id), non
    portee separement."""
    diff_weighted: AMSMeasurement   # code2="N-" : difference ponderee par `percent`, PAS normalisee
    diff_normalized: AMSMeasurement  # code2="S-" : 1 + difference ponderee, renormalisee (trace=3)
    diff_inverse: AMSMeasurement     # code2="SI" : tenseur 2 EXPRIME dans le repere du tenseur 1 (via son inverse)


def subtract_tensors(m1: AMSMeasurement, m2: AMSMeasurement, percent: float = 100.0) -> SubtractionResult:
    """Equivalent de `soustract` pour exactement 2 mesures (Util_AMS.f95:176-283) :
    `percent` (0-100, comme le prompt Fortran "% de soustraction entre 1 et
    100") pondere la contribution du 1er tenseur. Les 3 resultats reprennent
    l'id/etape/orientation du 2e tenseur (ams(2)), sauf `diff_weighted` qui
    prend le code2 fixe 'N-', `diff_normalized` 'S-', `diff_inverse` 'SI' -
    transcrit tel quel, y compris le fait que `diff_weighted` n'est PAS
    renormalise (contrairement aux deux autres) - verifie contre le source,
    pas une simplification."""
    pct = percent / 100.0

    n_minus = AMSMeasurement(
        id=m2.id, etape=m2.etape, code2="N-", cin=m2.cin, caz=m2.caz,
        dip=m2.dip, str_=m2.str_, vol=m2.vol,
        k11=(m1.s * m1.k11) - pct * (m2.s * m2.k11),
        k22=(m1.s * m1.k22) - pct * (m2.s * m2.k22),
        k33=(m1.s * m1.k33) - pct * (m2.s * m2.k33),
        k12=m1.s * m1.k12 - pct * m2.s * m2.k12,
        k23=m1.s * m1.k23 - pct * m2.s * m2.k23,
        k13=m1.s * m1.k13 - pct * m2.s * m2.k13,
    )
    n_minus.s = (n_minus.k11 + n_minus.k22 + n_minus.k33) / 3.0

    s_minus = AMSMeasurement(**{**m2.__dict__})
    s_minus.code2 = "S-"
    s_minus.k11 = 1 + m2.k11 - pct * m1.k11
    s_minus.k22 = 1 + m2.k22 - pct * m1.k22
    s_minus.k33 = 1 + m2.k33 - pct * m1.k33
    s_minus.k12 = m2.k12 - pct * m1.k12
    s_minus.k23 = m2.k23 - pct * m1.k23
    s_minus.k13 = m2.k13 - pct * m1.k13
    rnorm = (s_minus.k11 + s_minus.k22 + s_minus.k33) / 3.0
    s_minus.k11 /= rnorm
    s_minus.k12 /= rnorm
    s_minus.k13 /= rnorm
    s_minus.k22 /= rnorm
    s_minus.k23 /= rnorm
    s_minus.k33 /= rnorm
    s_minus.s = m2.s

    a, b, c, d, e, f = m1.k11, m1.k12, m1.k13, m1.k22, m1.k23, m1.k33
    det = a * d * f - a * e * e - b * b * f + 2 * b * c * e - c * c * d
    a11, a12, a13 = (d * f - e * e) / det, -(b * f - c * e) / det, (b * e - c * d) / det
    a21, a22, a23 = -(b * f - c * e) / det, (a * f - c * c) / det, -(a * e - b * c) / det
    a31, a32, a33 = (b * e - c * d) / det, -(a * e - b * c) / det, (a * d - b * b) / det
    s_i = AMSMeasurement(**{**m2.__dict__})
    s_i.code2 = "SI"
    s_i.k11 = m2.k11 * a11 + m2.k12 * a21 + m2.k13 * a31
    s_i.k12 = m2.k11 * a12 + m2.k12 * a22 + m2.k13 * a32
    s_i.k13 = m2.k11 * a13 + m2.k12 * a23 + m2.k13 * a33
    s_i.k22 = m2.k22 * a22 + m2.k23 * a32
    s_i.k23 = m2.k22 * a23 + m2.k23 * a33
    s_i.k33 = m2.k33 * a33

    return SubtractionResult(diff_weighted=n_minus, diff_normalized=s_minus, diff_inverse=s_i)
