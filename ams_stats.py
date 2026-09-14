"""
Statistiques tensorielles AMS : port de `anisotropie.f` (elprop/param/erbar/
ellips/tsmean/mds/fisnec), d'apres le rapport d'exploration verifie contre de
vraies donnees (magicams.txt/magicams_results.txt, site "98PL04" - formules
L/F/P/P' confirmees exactes, T confirme a la precision d'affichage pres).

Choix delibere : `elprop` (decomposition propre analytique par methode
trigonometrique de Cardano) est remplace par `numpy.linalg.eigh` (meme
resultat physique pour une matrice symetrique reelle, deja le choix retenu
cote STARpaleomag_Py pour l'ACP - voir calcul.linear_fit) plutot que de retranscrire
le solveur cubique du Fortran. En revanche, les SEUILS de classification
isotrope/prolate/oblate/triaxial (utilises par `tsmean` pour exclure les
tenseurs isotropes de la moyenne) sont repris EXACTEMENT du Fortran, y
compris leur asymetrie (0.00001/0.000001/0.00001 - pas une coquille a
harmoniser, confirme par lecture directe du source)."""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from ams_selection import AMSMeasurement, apply_orientation

# Hext/Fisher F(2,n-2) 95% table (anisotropie.f, fonction `fisnec`, 40
# entrees verbatim) - indexe 1..40 cote Fortran (fisnec(n-2)), ici stocke
# 0-indexed (FISNEC[0] == fisnec(1) du Fortran).
FISNEC = [
    199.5, 19.0, 9.55, 6.94, 5.79, 5.14, 4.74, 4.46, 4.26,
    4.1, 3.98, 3.89, 3.81, 3.74, 3.68, 3.63, 3.59, 3.55, 3.52, 3.49, 3.47,
    3.44, 3.42, 3.4, 3.39, 3.37, 3.35, 3.34, 3.33, 3.32, 3.31, 3.3, 3.29, 3.28,
    3.27, 3.26, 3.25, 3.24, 3.23, 3.22,
]


def fisnec(n: int) -> float:
    """Equivalent de `fisnec(i)` avec `i=n-2` (convention d'appel de
    `ellips`) : F(2,n-2) a 95% - `n>42` -> approximation asymptotique 3.0
    (meme seuil que le Fortran)."""
    if n > 42:
        return 3.0
    idx = n - 2
    if idx < 1:
        idx = 1
    if idx > 40:
        idx = 40
    return FISNEC[idx - 1]


def tensor_matrix(k11: float, k22: float, k33: float, k12: float, k23: float, k13: float) -> np.ndarray:
    return np.array([[k11, k12, k13], [k12, k22, k23], [k13, k23, k33]], dtype=float)


def eigen_decompose(a: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Equivalent de `elprop` (anisotropie.f:1097-1208) pour la partie
    calcul (PAS le solveur cubique lui-meme, remplace par `numpy.linalg.eigh`
    - voir docstring module). Retourne (valeurs propres decroissantes
    av1>=av2>=av3, vecteurs propres en COLONNES dans le meme ordre) - meme
    convention que le Fortran ou `elprop` garantit cet ordre par construction
    de la resolution trigonometrique (confirme par le rapport d'exploration,
    pas suppose)."""
    w, v = np.linalg.eigh(a)  # ascendant
    order = np.argsort(w)[::-1]
    return w[order], v[:, order]


def negative_k_corrected_tensor(a: np.ndarray) -> np.ndarray:
    """Corrige un tenseur ORIENTE `a` (voir `ams_selection.apply_
    orientation`) d'une mesure a susceptibilite negative (`m.s < 0` -
    tenseur entier deja retourne par Agico, code2="I-") AVANT de le
    moyenner avec d'autres (voir app.ouvrir_tsmean_dialog) : la direction
    associee a la plus grande valeur propre d'un tel tenseur correspond en
    pratique a kmin_phys (et inversement) - meme constat empirique que
    `ams_stats.principal_axes` (voir sa docstring, verifie sur 15 paires
    Re/Im reelles de 24WH.asc), demande explicite utilisateur ("les
    tenseurs moyens de Im negatifs ont aussi le kmax et le kmin
    intervertis").

    Pour un AFFICHAGE (`principal_axes`), simplement RELABELISER les
    directions apres decomposition suffit - mais pour une MOYENNE, cette
    correction doit se faire sur le TENSEUR lui-meme, AVANT de sommer/
    moyenner (demande explicite utilisateur - "comment faire une moyenne
    pour un site avec une partie des echantillons avec des Kim negatifs
    et d'autres positifs") : un site melangeant des specimens negatifs et
    positifs a un axe "grande valeur propre" qui ne represente PAS le
    meme role physique (kmax pour les positifs, kmin pour les negatifs)
    d'un specimen a l'autre - moyenner les tenseurs bruts tels quels
    reviendrait a additionner des axes incoherents entre eux. Reconstruit
    le tenseur avec les MEMES vecteurs propres mais les valeurs propres
    extremes echangees (`av[0]<->av[2]`, `av[1]`/kint inchange - meme
    invariant que `principal_axes`) : `V @ diag(av3,av2,av1) @ V.T`.
    Pour un groupe HOMOGENE (tout negatif, ou tout positif) le resultat
    est equivalent (aux imprecisions numeriques pres) a moyenner sans
    corriger puis ne relabeliser que la moyenne finale - mais seule cette
    version, appliquee mesure par mesure, reste correcte pour un groupe
    MELANGE."""
    av, vecs = eigen_decompose(a)
    av2 = av[[2, 1, 0]]
    return vecs @ np.diag(av2) @ vecs.T


def classify_ellipsoid(av: np.ndarray) -> int:
    """Equivalent de la classification isotrope/prolate/oblate/triaxial de
    `elprop` (seuils EXACTS, y compris leur asymetrie - anisotropie.f:1130-1150) :
    -1 triaxial, 0 isotrope, 1 prolate (revolution, k2~k3), 3 oblate
    (revolution, k1~k2). `av` : valeurs propres decroissantes, PEUVENT etre
    non normalisees (la classification normalise elle-meme par trace/3, comme
    le fait le Fortran avant ce test)."""
    rr = float(np.sum(av)) / 3.0
    if rr == 0:
        return 0
    a1, a2, a3 = av[0] / rr, av[1] / rr, av[2] / rr
    dif1 = (a1 - a3) / (a1 + a3) if (a1 + a3) else 0.0
    if dif1 < 0.00001:
        return 0
    dif2 = (a1 - a2) / (a1 + a2) if (a1 + a2) else 0.0
    if dif2 < 0.000001 or a2 < a1 / 200.0:
        return 3
    dif3 = (a2 - a3) / (a2 + a3) if (a2 + a3) else 0.0
    if dif3 < 0.00001:
        return 1
    return -1


def shape_params(k1: float, k2: float, k3: float) -> dict:
    """Equivalent de `param`/`paramres`/`paramech` (anisotropie.f:1804-1877,
    3336-3380) : k1>=k2>=k3 valeurs propres (unites physiques). Formules
    verifiees EXACTES contre magicams.txt reel (P/P'/L/F, site 98PL0401A) :
    P'=exp(sqrt(2*sum((ln(ki)-mean)^2))), T=(2*eta2-eta1-eta3)/(eta1-eta3).

    k1/k2/k3 non tous strictement positifs (susceptibilite negative -
    voir ams_asc.py, code2="I-" : le tenseur entier est retourne mais un
    axe individuel peut rester negatif, confirme sur un vrai specimen
    KLY5 de 24WH.asc) rendent alog() indefini - en Fortran comme en
    Python (`math.log` d'un negatif : ValueError, crash reel reproduit
    sur "List measurements"/Flinn/XY plots des que la selection contient
    un tel specimen). Retourne NaN partout plutot que de planter -
    l'appelant affiche alors "nan"/"n/a" pour ce point au lieu de perdre
    toute la liste/le graphe."""
    if k1 <= 0.0 or k2 <= 0.0 or k3 <= 0.0:
        nan = float("nan")
        return {
            "k0": nan, "km": (k1 + k2 + k3) / 3.0, "pan": nan, "P": nan, "Pprim": nan, "T": nan,
            "L": nan, "F": nan, "am1": nan, "am2": nan, "am3": nan, "r": nan, "ak": nan,
        }
    k0 = (k1 * k2 * k3) ** (1.0 / 3.0)
    km = (k1 + k2 + k3) / 3.0
    pan = (k1 - k3) / k3 * 100.0
    u1, u2, u3 = math.log(k1), math.log(k2), math.log(k3)
    u = (u1 + u2 + u3) / 3.0
    pprim = math.exp(math.sqrt(2.0 * ((u1 - u) ** 2 + (u2 - u) ** 2 + (u3 - u) ** 2)))
    t = (2 * u2 - u1 - u3) / (u1 - u3) if (u1 - u3) else 0.0
    am1, am2, am3 = (k1 - k0) / k0, (k2 - k0) / k0, (k3 - k0) / k0
    r = k1 / k2 + k2 / k3 - 1.0
    # ak (facteur de forme Flinn/Ramsay) : seuil et sentinelle EXACTS du
    # Fortran (anisotropie.f:1817-1822 / 3355-3360) - pas une approximation.
    if (k2 / k3 - 1.0) > 0.00001:
        ak = (k1 / k2 - 1.0) / (k2 / k3 - 1.0)
    else:
        ak = 1.0e5
    return {
        "k0": k0, "km": km, "pan": pan, "P": k1 / k3, "Pprim": pprim, "T": t,
        "L": k1 / k2, "F": k2 / k3,
        "am1": am1, "am2": am2, "am3": am3, "r": r, "ak": ak,
    }


def _rpc(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Equivalent de `rpc` (anisotropie.f:1241-1248)."""
    dec = math.degrees(math.atan2(y, x))
    if dec < 0.0:
        dec += 360.0
    inc = math.degrees(math.atan2(z, math.hypot(x, y)))
    ro = math.sqrt(x * x + y * y + z * z)
    return ro, dec, inc


@dataclass
class AxisResult:
    """Un axe (k1, k2 ou k3) d'un tenseur moyen : valeur propre, ecart-type
    (Jelinek), declinaison/inclinaison, et l'ellipse de confiance a 95%
    (matrice 2x2 `restm`, angles `alpha` = (demi-grand axe, demi-petit axe)
    en degres) - equivalent d'une "colonne" de `amsresults.axk1`/`tmres`/
    `alphares`."""
    eigenvalue: float = 0.0
    std: float = 0.0
    dec: float = 0.0
    inc: float = 0.0
    restm: Optional[np.ndarray] = None  # 2x2
    alpha: Tuple[float, float] = (0.0, 0.0)  # (alpha1=grand axe, alpha2=petit axe)


@dataclass
class TensorialMeanResult:
    """Equivalent de `amsresults` (AMS_OSX.inc) pour UN tenseur moyen : les
    3 axes (k1>=k2>=k3), le nombre de specimens effectivement utilises
    (`n`, apres exclusion des tenseurs isotropes), et le tenseur moyen brut
    (utile pour l'export/la correction d'anisotropie). `id` : PAS peuple
    par `tsmean` elle-meme (qui ne recoit que des tenseurs bruts, pas les
    mesures) - a renseigner par l'appelant, meme convention que `storeres`
    (anisotropie.f:2989-3049) : 6 premiers caracteres du 1er specimen +
    6 premiers du dernier."""
    n: int = 0
    axes: List[AxisResult] = field(default_factory=list)  # [k1, k2, k3]
    mean_tensor_normalized: Optional[np.ndarray] = None  # trace = 3
    ellipsoid_type: int = -1
    id: str = ""
    # "RE"/"IM" quand TOUTES les mesures moyennees venaient du meme canal
    # (voir app.ouvrir_tsmean_dialog, qui separe deja Re/Im - "pour le
    # calcul des tenseurs moyens, il vaut mieux separer Re et Im") ; None
    # sinon (groupe mixte ou source non classifiee). Utilise par app.
    # _save_mean_results_to_pmagani pour FORCER le code2 ecrit dans le
    # .pmagani a "RE"/"IM" plutot que de le redemander a l'utilisateur -
    # demande explicite utilisateur ("force the code to Re and Im when
    # saving the results in the file").
    source_code2: Optional[str] = None
    # "Y"/"N" - inclus dans l'export MagIC ou non, colonne .pmagani
    # DEDIEE (meme convention/schema que AMSMeasurement.export - le
    # niveau specimen) - demande explicite utilisateur ("is it possible
    # to export from AMS_py only the data and mean tensors that we want
    # to export"). "Y" par defaut : un fichier jamais touche par
    # app.ouvrir_marquer_export_dialog continue d'exporter TOUTES ses
    # moyennes, comme avant cette fonctionnalite.
    export: str = "Y"


def erbar(eigenvectors: np.ndarray, cov6: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Equivalent de `erbar` (anisotropie.f:1605-1665) : propage la
    covariance 6x6 des composantes du tenseur (repere cartesien d'origine)
    dans le repere des vecteurs propres du tenseur moyen. `eigenvectors` :
    matrice 3x3, colonnes = vecteurs propres (ordre k1,k2,k3). `cov6` :
    covariance 6x6 des composantes (k11,k22,k33,k12,k23,k13), POPULATION
    (denominateur N, pas N-1 - meme convention que `tsmean`).

    Retourne (vp: covariance 6x6 dans le repere propre, eigen_std: les 3
    ecarts-types des valeurs propres = sqrt(vp[i,i]) pour i=0,1,2)."""
    p = eigenvectors
    t = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            k = (j + 1) % 3
            l = (i + 1) % 3
            t[i, j] = p[j, i] * p[j, i]
            t[i, j + 3] = p[j, i] * p[k, i] * 2.0
            t[j + 3, i] = p[i, j] * p[i, k]
            t[i + 3, j + 3] = p[j, i] * p[k, l] + p[k, i] * p[j, l]
    vp = t @ cov6 @ t.T
    eigen_std = np.array([math.sqrt(max(0.0, vp[i, i])) for i in range(3)])
    return vp, eigen_std


def ellips(av: np.ndarray, vp6: np.ndarray, n: int) -> Tuple[List[np.ndarray], List[Tuple[float, float]]]:
    """Equivalent de `ellips` (anisotropie.f:1516-1604) : ellipse de
    confiance a 95% (Jelinek 1978) autour de chaque axe, a partir des
    valeurs propres `av` (decroissantes) et de la covariance en repere
    propre `vp6` (sortie de `erbar`). `n` : nombre de specimens (pour la
    table F de Hext). Retourne (restm[i]: matrice 2x2, alpha[i]:
    (alpha_grand_axe, alpha_petit_axe) en degres) pour i=0(k1),1(k2),2(k3)."""
    restm_list: List[np.ndarray] = []
    alpha_list: List[Tuple[float, float]] = []
    eps = 0.00001
    f = 3.0 if n > 42 else fisnec(n)
    c = 2.0 * (n - 1.0) * f / (n * (n - 2.0)) if n > 2 else 0.0

    w2x2 = []
    lambdas = []
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        denom_ij = (av[i] - av[j])
        denom_ik = (av[i] - av[k])
        a = vp6[3 + i, 3 + i] / (denom_ij * denom_ij) if denom_ij else 0.0
        b = vp6[3 + i, 3 + k] / (denom_ij * denom_ik) if (denom_ij and denom_ik) else 0.0
        c33 = vp6[3 + k, 3 + k] / (denom_ik * denom_ik) if denom_ik else 0.0
        w2x2.append(np.array([[a, b], [b, c33]]))
        delta = (a + c33) ** 2 - 4.0 * (a * c33 - b * b)
        delta = max(0.0, delta)
        l1 = ((a + c33) + math.sqrt(delta)) / 2.0
        l2 = ((a + c33) - math.sqrt(delta)) / 2.0
        lambdas.append((l1, l2))

    for i in range(3):
        a, b = w2x2[i][0, 0], w2x2[i][0, 1]
        cc = w2x2[i][1, 1]
        v = np.zeros((2, 2))
        for jx, lam in enumerate(lambdas[i]):
            lam_c = max(lam, 0.000001)
            dif1 = a - lam
            dif2 = cc - lam
            if abs(dif1) < eps:
                if abs(dif2) < eps:
                    v[:, jx] = [1.0, 0.0] if jx == 0 else [0.0, 1.0]
                    continue
                v1 = 1.0
                v2 = -dif1 / b if b else 0.0
            else:
                v2 = 1.0
                v1 = -dif2 / b if b else 0.0
            r = math.hypot(v1, v2)
            if r == 0:
                r = 1.0
            v[0, jx], v[1, jx] = v1 / r, v2 / r
        e = np.diag([math.sqrt(c * max(lambdas[i][0], 0.000001)),
                     math.sqrt(c * max(lambdas[i][1], 0.000001))])
        restm = v @ e @ v.T
        restm_list.append(restm)
        alpha1 = math.degrees(math.atan(math.sqrt(c * max(lambdas[i][0], 0.000001))))
        alpha2 = math.degrees(math.atan(math.sqrt(c * max(lambdas[i][1], 0.000001))))
        alpha_list.append((alpha1, alpha2))

    return restm_list, alpha_list


def tsmean(
    tensors: List[Tuple[float, float, float, float, float, float]],
    normalize: bool = True,
) -> Optional[TensorialMeanResult]:
    """Equivalent de `tsmean` (anisotropie.f:1668-1801) : moyenne
    tensorielle avec statistiques de Jelinek (1978). `tensors` : liste de
    (k11,k22,k33,k12,k23,k13) - normalises ou bruts selon `normalize`
    (True = chaque tenseur est divise par sa propre trace/3 avant d'etre
    moyenne, meme convention par defaut que le prompt Fortran
    "normalisation of the tensors (Y/n)?", defaut Y). Necessite au moins 3
    tenseurs (`il>2`, comme le Fortran - sinon renvoie None, equivalent du
    message "nombre de donnees insuffisant"). Les tenseurs isotropes (selon
    `classify_ellipsoid`) sont EXCLUS de la moyenne, comme le Fortran."""
    if len(tensors) <= 2:
        return None

    used = []
    for k in tensors:
        a = tensor_matrix(*k)
        av, _ = eigen_decompose(a)
        if classify_ellipsoid(av) == 0:
            continue  # tenseur isotrope exclu, comme le Fortran
        used.append(k)
    n = len(used)
    if n < 1:
        return None

    sk = np.zeros(6)
    for k in used:
        r = (k[0] + k[1] + k[2]) / 3.0 if normalize else 1.0
        sk += np.array(k) / r
    sk /= n

    cov6 = np.zeros((6, 6))
    for k in used:
        r = (k[0] + k[1] + k[2]) / 3.0 if normalize else 1.0
        d = np.array(k) / r - sk
        cov6 += np.outer(d, d)
    cov6 /= n  # population (denominateur N), pas N-1 - confirme par le rapport

    a_mean = tensor_matrix(*sk)
    av, vecs = eigen_decompose(a_mean)
    ellipsoid_type = classify_ellipsoid(av)

    result = TensorialMeanResult(n=n, mean_tensor_normalized=a_mean, ellipsoid_type=ellipsoid_type)
    if ellipsoid_type == 0:
        # tenseur moyen isotrope : pas de statistiques d'axes (le Fortran
        # retourne directement dans ce cas, "601 write(*,999); return").
        return result

    vp6, eigen_std = erbar(vecs, cov6)
    restm_list, alpha_list = ellips(av, vp6, n)

    for i in range(3):
        x, y, z = vecs[0, i], vecs[1, i], vecs[2, i]
        _, dec, inc = _rpc(x, y, z)
        result.axes.append(AxisResult(
            eigenvalue=float(av[i]), std=float(eigen_std[i]), dec=dec, inc=inc,
            restm=restm_list[i], alpha=alpha_list[i],
        ))
    return result


def mean_susceptibility(bulk_values: List[float]) -> Optional[dict]:
    """Equivalent de `mds` (anisotropie.f:3560-3623) : moyenne arithmetique
    et geometrique (+ ecart-type geometrique, log base 10) de valeurs de
    susceptibilite bulk. Necessite au moins 3 valeurs (meme garde que le
    Fortran, `if(itot==1) return` / `if(itot==2) return`)."""
    n = len(bulk_values)
    if n < 3:
        return None
    arith_mean = sum(bulk_values) / n
    logs = [math.log10(v) for v in bulk_values]
    rat1 = sum(logs) / n
    somcar = sum((lv - rat1) ** 2 for lv in logs)
    ecart = math.sqrt(somcar / (n - 1))  # N-1, PAS N (different de tsmean - confirme)
    return {
        "arith_mean": arith_mean,
        "geom_mean": 10 ** rat1,
        "geom_sup": 10 ** (rat1 + ecart),
        "geom_inf": 10 ** (rat1 - ecart),
        "n": n,
    }


def principal_axes(
    m: AMSMeasurement, orientation: int, flip_negative_inclination: bool = True,
) -> List[Tuple[float, float, float]]:
    """Equivalent du peuplement de `axes(j,1,i)/axes(j,2,i)/axes(j,3,i)` pour
    UNE mesure (`editn`/`editams`, anisotropie.f:732-859/2747-2869) : decompose
    le tenseur ORIENTE de `m` (`selection.apply_orientation`) et retourne les
    3 axes k1>=k2>=k3 sous la forme (valeur_propre, declinaison,
    inclinaison).

    `flip_negative_inclination` (defaut True) : ramene l'inclinaison
    toujours positive (retournement antipodal dec+180 mod 360) - meme
    normalisation que `liste`/`lisresmem`, appliquee ICI a chaque appel
    (fonction pure) plutot que mutee une fois pour toutes dans un tableau
    partage comme le fait le Fortran. Mettre a False pour le stereonet
    (`ams_stereo`), qui a SA PROPRE logique de retournement (`iinv`,
    distincte de celle de `liste` - un axe d'inclinaison negative y garde
    sa position reelle avec un symbole different si `iinv` n'est pas
    actif, transcrit de `stereo`/`stereo2`).

    `m.s < 0` (susceptibilite negative - voir ams_asc.py, tenseur entier
    deja retourne par Agico, code2="I-") : le vecteur associe a la plus
    grande valeur propre (apres ce retournement) correspond en pratique a
    la direction kmin du signal REEL du meme specimen, et inversement -
    verifie directement sur 15 paires Re/Im reelles de 24WH.asc (chaque
    axe Im s'aligne avec l'axe Re oppose une fois les DIRECTIONS k1/k3
    echangees) - demande explicite utilisateur ("flip the vector kmax
    and kmin for negative susceptibility ... more coherent with what is
    observed for the real part"). Seules les DIRECTIONS (dec/inc) sont
    echangees entre les positions 0 et 2 - les valeurs propres restent a
    leur position (tri decroissant inchange), pour que L/F/P/T/P'
    (ams_stats.shape_params, purs ratios de valeurs propres) restent
    calcules exactement comme avant ; seul le libelle "quel axe est
    kmax/kmin" change pour ce cas."""
    a = apply_orientation(m, orientation)
    av, vecs = eigen_decompose(a)
    dirs = []
    for i in range(3):
        x, y, z = vecs[0, i], vecs[1, i], vecs[2, i]
        _, dec, inc = _rpc(x, y, z)
        if flip_negative_inclination and inc < 0.0:
            inc = -inc
            dec += 180.0
            if dec > 360.0:
                dec -= 360.0
        dirs.append((dec, inc))
    if m.s < 0.0:
        dirs[0], dirs[2] = dirs[2], dirs[0]
    return [(float(av[i]), dirs[i][0], dirs[i][1]) for i in range(3)]


def format_measurement_list(measurements: List[AMSMeasurement], orientation: int) -> str:
    """Equivalent de `liste`/`lismes` (anisotropie.f:1250-1403, format 1404) :
    PAS le tenseur brut du fichier .ANI (ca, c'est `listeani`/"Liste File
    .ANI") mais les axes propres k1/k2/k3 (dec/inc) + parametres de forme
    L/F/P/P'/T, calcules a partir du tenseur ORIENTE (selon l'orientation
    couramment choisie) de chaque mesure. Susceptibilite affichee : `s*1e-5`
    sauf si le code2 commence par 'A' ou 'B', auquel cas c'est l'axe k2 -
    quirk reel du Fortran, transcrit tel quel (pas une approximation).

    Regroupe les lignes par specimen (comme `sterams2`, la vue "un seul
    specimen a la fois" du menu Stereographic sample - meme regroupement
    par id consecutif) et affiche, juste APRES le groupe de lignes d'un
    specimen, le commentaire `info` de sa ligne 'A0' UNE SEULE fois (pas
    repete sur chacune des 15 variantes jackknife, qui partagent toutes le
    meme texte de base) - demande explicite utilisateur ("can you write
    the info comment from the A0 line in the text window just below the
    listing of the anisotropy parameters"). 'A0' est la ligne qui porte
    le verdict PmagPy satisfactory/not satisfactory (voir calcul.
    _format_pmagani_line cote STARpaleomag_Py) ; les autres variantes n'en ont
    pas d'equivalent propre."""
    lines = []
    group: List[AMSMeasurement] = []

    def _flush_group():
        if not group:
            return
        for m in group:
            axes = principal_axes(m, orientation)
            k1, k2, k3 = axes[0][0], axes[1][0], axes[2][0]
            sp = shape_params(k1, k2, k3)
            code1 = (m.code2 or "  ")[0:1]
            code2c = (m.code2 or "  ")[1:2]
            valsuscnrm = k2 if code1 in ("A", "B") else m.s * 1.0e-5
            lines.append(
                f"{m.id:<12} {m.etape:4d} {code1}{code2c} {valsuscnrm:12.3E}  "
                f"{sp['L']:6.3f}  {sp['F']:6.3f}  {sp['P']:6.3f}    "
                f"{axes[0][1]:5.1f} {axes[0][2]:5.1f}    {axes[1][1]:5.1f} {axes[1][2]:5.1f}    "
                f"{axes[2][1]:5.1f} {axes[2][2]:5.1f}   {sp['T']:7.3f}  {sp['Pprim']:6.3f}"
            )
        a0 = next((m for m in group if (m.code2 or "").strip() == "A0" and m.info), None)
        if a0 is not None:
            # coupe en deux lignes au marqueur PmagPy (voir calcul.
            # _format_pmagani_line cote STARpaleomag_Py, meme separateur " - ")
            # - demande explicite utilisateur ("split on two lines") : la
            # ligne combinee est trop longue pour rester lisible en un
            # seul morceau.
            marker = " - PmagPy Hext F-test:"
            if marker in a0.info:
                base, _sep, verdict = a0.info.partition(marker)
                lines.append(f"  {base.strip()}")
                lines.append(f"  PmagPy Hext F-test:{verdict}")
            else:
                lines.append(f"  {a0.info}")

    for m in measurements:
        if group and group[-1].id != m.id:
            _flush_group()
            group = []
        group.append(m)
    _flush_group()

    return "\n".join(lines) + ("\n" if lines else "")


def fortran_e(x: float, digits: int = 7) -> str:
    """Format Fortran classique E14.7 : mantisse normalisee 0.dddddddE±dd
    (PAS la notation scientifique Python usuelle, ex. 1.0095790 -> Fortran
    "0.1009579E+01" contre Python "1.0095790E+00") - necessaire pour
    reproduire l'affichage reel de `tsmean`."""
    if x == 0:
        return "0." + "0" * digits + "E+00"
    sign = "-" if x < 0 else ""
    x = abs(x)
    exp = math.floor(math.log10(x)) + 1
    mantissa = round(x / (10.0 ** exp), digits)
    if mantissa >= 1.0:
        mantissa /= 10.0
        exp += 1
    elif mantissa < 0.1:
        mantissa *= 10.0
        exp -= 1
    mant_str = f"{mantissa:.{digits}f}"
    exp_sign = "+" if exp >= 0 else "-"
    return f"{sign}{mant_str}E{exp_sign}{abs(exp):02d}"


def _box_border() -> str:
    return " " + "*" * 63


def _box_line(content: str = "") -> str:
    return " *" + content[:61].ljust(61) + "*"


_ORIENTATION_BOX_LABEL = {
    # "geographic coordinates" (iorient=2) est confirme sur un vrai
    # transcript utilisateur reel (tensorial mean, 11 mesures). Les
    # libelles pour 1/3 sont deduits par analogie (non retrouves verbatim
    # dans le source ni verifies sur un transcript reel) - a corriger si
    # un vrai cas les infirme.
    1: "sample coordinates",
    2: "geographic coordinates",
    3: "tectonic coordinates",
}


def format_tsmean_box(result: TensorialMeanResult, orientation: int = 2) -> str:
    """Reproduit la boite decorative reelle affichee par `tsmean`/`param`
    apres calcul de la moyenne tensorielle - largeurs de colonnes et
    formats mesures directement sur un transcript utilisateur reel (11
    mesures, coordonnees geographiques), verifies caractere pres (22/23
    lignes identiques, le seul ecart residuel etant un artefact
    d'arrondi lie aux valeurs d'entree de test, pas un defaut de
    formatage)."""
    if result.ellipsoid_type == 0 or not result.axes:
        return " Isotropic mean tensor - no axis statistics.\n"

    k1, k2, k3 = (ax.eigenvalue for ax in result.axes)
    bulk = (k1 + k2 + k3) / 3.0
    sp = shape_params(k1, k2, k3)
    label = _ORIENTATION_BOX_LABEL.get(orientation, "geographic coordinates")

    lines = [" Parameters of the tensorial mean :", "", _box_border(), _box_line()]
    lines.append(_box_line("  " + label.ljust(34) + "bulk= " + fortran_e(bulk)))
    lines.append(_box_line())
    for i, ax in enumerate(result.axes, start=1):
        prefix = f"  k{i}= {fortran_e(ax.eigenvalue)}".ljust(24)
        mid = (prefix + f"dec={ax.dec:5.1f}").ljust(36)
        lines.append(_box_line(mid + f"inc={ax.inc:5.1f}"))
    lines.append(_box_line())
    lines.append(_box_border())
    lines.append(_box_line(" number of data =" + f"{result.n:5d}"))
    lines.append(_box_line())
    lines.append(_box_line(
        "   lineation l: k1/k2=" + f"{sp['L']:7.3f}"
        + "     mi=(ki-ko)/ko:  " + f"m1={sp['am1']:6.3f}  "))
    lines.append(_box_line(
        "   foliation f: k2/k3=" + f"{sp['F']:7.3f}" + " " * 21 + f"m2={sp['am2']:6.3f}  "))
    lines.append(_box_line(
        "   degree    p: k1/k3=" + f"{sp['P']:7.3f}" + " " * 21 + f"m3={sp['am3']:6.3f}  "))
    lines.append(_box_line())
    lines.append(_box_border())
    lines.append(_box_line())
    shape_label = "oblate" if sp["ak"] <= 1.0 else "prolate"
    lines.append(_box_line(
        "   shape:      k=" + f"{sp['ak']:7.3f}" + "     t=" + f"{sp['T']:5.2f}"
        + " " * 7 + f"{shape_label:<18}"))
    lines.append(_box_line(
        "   intensity:  r=" + f"{sp['r']:7.3f}" + "     %an=" + f"{sp['pan']:5.1f}"
        + "     p'=" + f"{sp['Pprim']:6.3f}" + " " * 9))
    lines.append(_box_line())
    lines.append(_box_border())
    return "\n".join(lines) + "\n"


_ORIENT_CODE = {1: "Sa", 2: "IS", 3: "TC"}  # meme convention que app._ORIENT_LABEL


def _lisresmem_row(index: int, orientation: int, result: TensorialMeanResult) -> str:
    """Une ligne du tableau `lisresmem` (voir format_lisresmem_table) -
    reproduit FORMAT 202 (anisotropie.f:3161-3162) : index, id, code
    orientation, n, puis pour chaque axe (k1=Kmax, k2=Kint, k3=Kmin)
    valeur propre/D/I/p1/p2, puis lineation/foliation/anisotropie totale,
    forme (k, t, oblate/prolate), intensite (r, Pprim)."""
    if result.ellipsoid_type == 0 or not result.axes:
        return f"{index:2d}: {result.id:<12s}  (isotropic tensor - no axis statistics)"

    k1, k2, k3 = (ax.eigenvalue for ax in result.axes)
    sp = shape_params(k1, k2, k3)
    shape_label = "oblate" if sp["ak"] <= 1.0 else "prolate"
    orient = _ORIENT_CODE.get(orientation, "??")

    row = f"{index:2d}: {result.id:<12s} {orient:<2s} {result.n:2d}"
    for ax in result.axes:
        p1, p2 = ax.alpha
        dec, inc = ax.dec, ax.inc
        if inc < 0.0:
            # display-only : Fortran (anisotropie.f:3141-3148) force
            # l'inclinaison positive et ajoute 180 a la declinaison avant
            # d'imprimer - PAS applique aux axes eux-memes ailleurs dans
            # le port (tsmean/_rpc ne le fait pas), seulement ici.
            inc = -inc
            dec = (dec + 180.0) % 360.0
        row += f"  {ax.eigenvalue:6.3f}  {dec:5.1f} {inc:5.1f} {p1:4.1f} {p2:4.1f}"
    row += (
        f"   {sp['L']:6.3f}  {sp['F']:6.3f}  {sp['P']:6.3f}"
        f" {sp['ak']:6.2f} {sp['T']:5.2f}  {shape_label:<8s}  {sp['r']:6.3f}  {sp['Pprim']:6.3f}"
    )
    return row


def format_lisresmem_table(entries) -> str:
    """Reproduit le tableau COMPACT (une ligne par site) de `lisresmem`
    (anisotropie.f:3094-3163, FORMAT 202) - demande explicite utilisateur
    ("is it possible to use the original form of listing the mean tensor
    per site in memory (check the Fortran source) and not the detailed
    one page information") : `format_tsmean_box` reproduit la boite
    decorative de `tsmean`/`param` (affichee juste APRES le calcul d'UN
    tenseur), pas celle de `lisresmem` - qui, source Fortran a l'appui,
    n'a jamais ete une boite mais un tableau, un site par ligne. `entries`
    : iterable de (orientation, TensorialMeanResult), meme forme que
    self.mean_results (app.py)."""
    header = (
        "                                 Kmax                         Kint"
        "                        kmin                                         shape"
        "               intensity\n"
        "                       ------------------------      "
        "------------------------      ------------------------                           "
        "---------------      ------------\n"
        " ID tenseur         N    k1     D1    I1   p1  p2      k2     D2    I2   p1"
        "  p2      k3     D3    I3   p1  p2      lin    fol    ani     k    t    forme"
        "     r     Pprim"
    )
    rows = [_lisresmem_row(i, orientation, res) for i, (orientation, res) in enumerate(entries, start=1)]
    return "\n".join([header] + rows) + "\n"


# orientation (1=Sa/2=IS/3=TC, voir ams_selection._ORIENT_MODE_TAG) -> code
# MagIC officiel de aniso_tilt_correction (sites.txt ET specimens.txt,
# meme colonne/memes codes) - VERIFIE contre le data model reel
# (data_model.json, tables "sites"/"specimens", colonne
# "aniso_tilt_correction" : "Correction between geographic (0%) and
# stratigraphic (100%); unoriented (-1%) ...") : -1 pour les coordonnees
# specimen, PAS "1" comme le code interne STARpaleomag_Py/calcul.
# _ORIENT_TO_FILE_CODE (colonne "IS/TC" de .pmagres) - deux conventions
# INDEPENDANTES, ne pas confondre. `aniso_tilt_correction` est
# `requiredIfGroup("Anisotropy")` dans le data model : une ligne
# anisotropie sans cette colonne echouerait la validation MagIC.
MAGIC_TILT_CORRECTION_CODE = {1: "-1", 2: "0", 3: "100"}


def magic_site_aniso_fields(
    result: TensorialMeanResult, aniso_type: str = "AMS", orientation: Optional[int] = None,
) -> Dict[str, str]:
    """Construit les colonnes MagIC `sites.txt` (data model 3, table
    "sites", groupe "Anisotropy") pour un tenseur moyen DEJA calcule
    (`TensorialMeanResult`, ex. par `tsmean` ou importe via
    `ams_selection.read_ani_mean_results`) - demande explicite
    utilisateur ("check in magic how to export a mean tensor at site
    level"). Champs distinctifs du niveau SITE (contrairement au niveau
    specimen, qui exporte `aniso_s`, voir app.exporter_magic_dialog) :
    `aniso_v1`/`aniso_v2`/`aniso_v3` (Matrix, un par axe propre - verifie
    dans pmagpy.ipmag, ligne ~10725 : "tau : dec : inc : eta/zeta : ...")
    plutot qu'un tenseur brut, et une dizaine de parametres de forme
    scalaires (aniso_p/l/f/t/pp/perc/total/ll/ff/vg/fl - tous verifies
    contre la description officielle de chaque colonne, MagIC-data-
    model.txt, groupe "sites"/"Anisotropy").

    `orientation` (1=Sa/2=IS/3=TC, meme convention que self.orientation/
    ams_selection._ORIENT_MODE_TAG) - demande explicite utilisateur
    ("are there other parameters that were not previously calculated and
    useful for the magic export") : ajoute `aniso_tilt_correction`
    (voir MAGIC_TILT_CORRECTION_CODE), colonne REQUISE par le data model
    des qu'un champ du groupe "Anisotropy" est renseigne
    (`requiredIfGroup("Anisotropy")`) et absente de cet export jusqu'ici.
    None (par defaut, ex. resultat sans orientation connue) omet la cle
    plutot que de deviner.

    eta/zeta (orientation + demi-angle de l'ellipse de confiance a 95%
    autour de chaque axe, vers les DEUX autres) : pmagpy lui-meme (dans
    ipmag.py, meme fonction que ci-dessus) pointe eta/zeta directement
    vers les AUTRES vecteurs propres (ex. l'eta de v1 = la direction de
    v2 elle-meme, PAS une orientation d'ellipse recalculee) plutot que
    de calculer une orientation d'ellipse pleinement libre - AMS_Py suit
    la MEME simplification ici, avec les demi-angles de Jelinek 1978
    deja calcules par `ellips()` (alpha1=demi-grand axe, alpha2=demi-
    petit axe - PAS les angles de paire e12/e13/e23 de Hext 1963 que
    pmagpy calcule separement pour un test a 6 positions specimen par
    specimen, un formalisme different, non recalcule ici) : ecart
    documente explicitement dans la colonne `description` de chaque
    ligne plutot que tu par une fausse equivalence.

    BUG REEL repere en verifiant cette convention (PAS reproduit ici) :
    dans pmagpy.ipmag (~ligne 10833), le bloc HextRec['aniso_v3'] utilise
    `hpars["e12"]` pour son propre angle "eta" (vers v1), alors que
    dohext() (pmag.py:8578) definit explicitement `e13` comme l'angle
    entre l'axe principal (v1) et l'axe mineur (v3) - v1 et v3 sont
    exactement les deux axes en jeu ici, `e13` est donc la valeur
    attendue, pas `e12` (visiblement un copier-coller depuis le bloc
    v1). Sans jeu de donnees reel pour confirmer si c'est une erreur de
    frappe ou un choix delibere non documente, AMS_Py applique ici la
    correspondance mathematiquement coherente avec la docstring de
    dohext plutot que de reproduire ce qui ressemble a une erreur.

    `aniso_type` : "AMS"/"AARM"/"ATRM"/"AIRM" (voir app._ANI_MAGIC_INFO,
    meme convention que l'export specimen). Retourne un dict pret a
    devenir une ligne de `sites.txt` (colonnes deja formatees en
    chaines) - specimen/sample/site/method_codes restent a la charge de
    l'appelant (pas connus de `TensorialMeanResult` lui-meme)."""
    k1, k2, k3 = (ax.eigenvalue for ax in result.axes)
    sp = shape_params(k1, k2, k3)

    # aniso_perc/aniso_total : formules officielles (MagIC-data-model.txt,
    # PAS deja calculees par shape_params - "pan" y est (k1-k3)/k3*100,
    # une definition VOISINE mais differente de aniso_perc/aniso_total).
    total_k = k1 + k2 + k3
    aniso_perc = 100.0 * (k1 - k3) / total_k if total_k else 0.0
    mean_k = total_k / 3.0
    aniso_total = 100.0 * (k1 - k3) / mean_k if mean_k else 0.0
    # aniso_ll/aniso_ff (Woodcock 1977, log-based - distinct de L/F de
    # Jelinek deja dans shape_params) et aniso_vg (Graham) : formules
    # officielles directes, pas encore ailleurs dans ce module.
    aniso_ll = math.log(k1 / k2) if k2 else 0.0
    aniso_ff = math.log(k2 / k3) if k3 else 0.0
    aniso_vg = math.degrees(math.asin(math.sqrt((k2 - k3) / (k1 - k3)))) if (k1 - k3) else 0.0
    aniso_fl = sp["F"] / sp["L"] if sp["L"] else 0.0

    dec = [ax.dec for ax in result.axes]
    inc = [ax.inc for ax in result.axes]
    alpha_major = [ax.alpha[0] for ax in result.axes]
    alpha_minor = [ax.alpha[1] for ax in result.axes]

    def aniso_v(i: int, j: int, k: int) -> str:
        # axe i, "eta" vers l'axe j (demi-grand axe de Jelinek), "zeta"
        # vers l'axe k (demi-petit axe) - voir docstring pour la
        # simplification (directions = les AUTRES vecteurs propres,
        # comme pmagpy, demi-angles = Jelinek alpha1/alpha2, pas Hext).
        return " : ".join(str(v) for v in [
            result.axes[i].eigenvalue, dec[i], inc[i], "eta/zeta",
            dec[j], inc[j], alpha_major[i],
            dec[k], inc[k], alpha_minor[i],
        ])

    fields = {
        "aniso_type": aniso_type,
        "aniso_v1": aniso_v(0, 1, 2),
        "aniso_v2": aniso_v(1, 0, 2),
        "aniso_v3": aniso_v(2, 0, 1),
        "aniso_p": f"{sp['P']:.6f}",
        "aniso_pp": f"{sp['Pprim']:.6f}",
        "aniso_t": f"{sp['T']:.6f}",
        "aniso_l": f"{sp['L']:.6f}",
        "aniso_f": f"{sp['F']:.6f}",
        "aniso_perc": f"{aniso_perc:.4f}",
        "aniso_total": f"{aniso_total:.4f}",
        "aniso_ll": f"{aniso_ll:.6f}",
        "aniso_ff": f"{aniso_ff:.6f}",
        "aniso_vg": f"{aniso_vg:.4f}",
        "aniso_fl": f"{aniso_fl:.6f}",
        "description": (
            "eta/zeta point toward the other two eigenvectors with Jelinek (1978) "
            "confidence semi-angles (alpha1/alpha2 from AMS_Py's tsmean), not Hext "
            "(1963) pairwise e12/e13/e23 - see magic_site_aniso_fields docstring"
        ),
    }
    if orientation in MAGIC_TILT_CORRECTION_CODE:
        fields["aniso_tilt_correction"] = MAGIC_TILT_CORRECTION_CODE[orientation]
    return fields
