"""
Port Python du modele de donnees AMS (Anisotropy of Magnetic Susceptibility) :
equivalent de la structure `/amsdata/` (AMS_OSX.inc) et des routines de lecture/
selection de `lect_asc.f`/`Util_AMS.f95` (AMS_OSX_AWE, reference/AMS_OSX_AWE).

Une mesure AMS = UNE position/orientation de mesure au kappabridge pour un
specimen donne (etape = position 1..15 selon le protocole a 15 positions
standard AGICO, ou 0/192 pour un tenseur deja moyenne par l'instrument selon
les fichiers reels observes) : id, orientation du carottier (cin/caz),
orientation du pendage (dip/str), volume, susceptibilite scalaire `s`, et le
tenseur de susceptibilite normalise (k11,k22,k33,k12,k23,k13, trace=3 - PAS
brut comme les tenseurs STARpaleomag_Py ANI).
"""

import math
import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np


@dataclass
class AMSMeasurement:
    """Equivalent d'un enregistrement `/amsdata/` (AMS_OSX.inc:2-10)."""
    id: str = ""
    etape: int = 0
    code2: str = "  "  # 2 caracteres (ex 'N0', 'A0', 'I-')
    cin: float = 0.0    # pendage (plunge) de l'axe du carottier
    caz: float = 0.0    # azimut de l'axe du carottier
    dip: float = 0.0    # pendage de la stratification/foliation de reference
    str_: float = 0.0   # direction (strike) de la stratification de reference
    vol: float = 1.0
    s: float = 0.0       # susceptibilite scalaire (moyenne), unites SI*1e5 (voir readasc)
    k11: float = 1.0
    k22: float = 1.0
    k33: float = 1.0
    k12: float = 0.0
    k23: float = 0.0
    k13: float = 0.0
    info: str = ""
    # Champs .pmagani (voir ams_selection._read_pmagani_file / calcul.
    # AniTensor cote STARpaleomag_Py, MEME schema de colonnes) - demande
    # explicite utilisateur ("What will be the most complete .pmagani
    # file"). None/absent pour les lignes venant d'un ancien .ANI (format
    # D, ne les a jamais eues) ou quand non calcule.
    n_positions: Optional[int] = None
    sigma: Optional[float] = None
    ftest: Optional[float] = None
    ftest12: Optional[float] = None
    ftest23: Optional[float] = None
    # 'g'/'b'/None - verdict PmagPy (Hext F-test), voir calcul.AniTensor.
    # quality cote STARpaleomag_Py - MEME colonne .pmagani ("quality"), pour que
    # les deux applications ecrivent/lisent exactement le meme schema
    # (sinon un fichier ecrit par l'une desaligne la lecture par colonnes
    # fixes de l'autre) - demande explicite utilisateur ("la correction
    # d'anisotropie ne se fait pas si le message satisfactory or not
    # n'est pas dans le fichier pmagani").
    quality: Optional[str] = None
    # "Y"/"N" - inclus dans l'export MagIC ou non, colonne .pmagani
    # DEDIEE (meme convention que calcul.AniTensor.export cote
    # STARpaleomag_Py, MEME fichier partage) - demande explicite
    # utilisateur ("is it possible to export from AMS_py only the data
    # and mean tensors that we want to export and only these selected
    # data will be taken into account in the main export from
    # Starpaleomag") : "Y" par defaut (tout exporte, comme avant cette
    # fonctionnalite - voir ouvrir_marquer_export_dialog, qui est la
    # SEULE facon de mettre "N" sur une ligne) pour qu'un fichier jamais
    # touche par cette fonctionnalite continue de tout exporter sans
    # rien changer a son comportement actuel.
    export: str = "Y"
    # Champ applique (A/m) durant la mesure - colonne .pmagani DEDIEE
    # (meme convention/fichier partage que calcul.AniTensor.field cote
    # STARpaleomag_Py) - demande explicite utilisateur ("in the AMS
    # measurements, there is also an important parameter, the field
    # used (in A/m)... The MFK2 instrument allow measurements at
    # different field values") : le kappabridge MFK2 peut mesurer un
    # meme lot de specimens a plusieurs champs differents (200, 425,
    # 700, 50, 5 A/m... vus sur un vrai fichier .asc) - perdu jusqu'ici,
    # jamais ecrit dans .pmagani bien que deja PARSE depuis le .asc
    # (voir ams_asc.parse_asc_file, marqueur F1/F3). None si non
    # disponible (ancien .ANI, ou branche "susc." du .asc - non
    # verifiee, voir ams_asc.py).
    field: Optional[float] = None


def is_imaginary_component(m: AMSMeasurement) -> bool:
    """True si `m` vient du canal imaginaire (KLY5 AC susceptibility -
    voir ams_asc.py, `component_note`) plutot que du canal reel/standard.
    Detecte via `info`, PAS `code2` seul : un tenseur negatif ecrase
    toujours `code2` en "I-" que le composant d'origine soit reel ou
    imaginaire (ams_asc.py:227-229 - dans les faits jamais observe pour
    un Re reel, mais le champ ne le garantit pas) ; `info` porte lui la
    distinction "real component"/"imaginary component" AVANT cet
    ecrasement, donc reste fiable dans tous les cas. Utilise par
    ams_stereo/ams_xy pour styliser les points Im differemment des Re
    (fill gris + stroke colore) - demande explicite utilisateur ("dans
    les plots ... mettre les tenseurs d'Im avec un fill en gris et le
    stroke color rouge vert bleu ... pour differencier les Re des Im")."""
    return "imaginary component" in (m.info or "")


def cart(r: float, d: float, ai: float) -> Tuple[float, float, float]:
    """Equivalent de `cart` (anisotropie.f:978-984) : (rayon, declinaison,
    inclinaison en degres) -> (x,y,z) cartesien."""
    w = math.pi / 180.0
    x = r * math.cos(d * w) * math.cos(ai * w)
    y = r * math.sin(d * w) * math.cos(ai * w)
    z = r * math.sin(ai * w)
    return x, y, z


def polere(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Equivalent de `polere` (anisotropie.f:1012-1034) : cartesien -> (rayon,
    declinaison, inclinaison en degres). Meme convention que
    `selection.polere` du port STARpaleomag_Py (verifiee identique formule a
    formule)."""
    horiz = math.hypot(x, y)
    mag = math.hypot(horiz, z)
    if mag == 0.0:
        return 0.0, 0.0, 0.0
    dec = math.degrees(math.atan2(y, x))
    if dec < 0.0:
        dec += 360.0
    inc = math.degrees(math.atan2(z, horiz))
    return mag, dec, inc


def _tensor_matrix(m: AMSMeasurement) -> np.ndarray:
    return np.array([
        [m.k11, m.k12, m.k13],
        [m.k12, m.k22, m.k23],
        [m.k13, m.k23, m.k33],
    ], dtype=float)


def corfor_tensor(a: np.ndarray, cin: float, caz: float) -> np.ndarray:
    """Equivalent de `corfor` pour un TENSEUR (anisotropie.f:862-887, PAS le
    `corfor` vectoriel de `selection.py` cote STARpaleomag_Py - meme nom, signature
    differente) : rotation du tenseur `a` (matrice 3x3 symetrique, repere
    echantillon) vers le repere in-situ, via l'axe du carottier (pendage
    `cin`, azimut `caz`). Construit la matrice de changement de base `tm`
    exactement comme le Fortran (3 vecteurs cart() a partir de d=caz, ai=cin)
    puis applique la transformation par similarite b_ij = sum_km tm_ik tm_jm
    a_km."""
    d, ai = caz, cin
    tm = np.zeros((3, 3))
    tm[:, 0] = cart(1.0, d - 90.0, -ai)
    tm[:, 1] = cart(1.0, d, 0.0)
    tm[:, 2] = cart(1.0, d - 90.0, 90.0 - ai)
    return tm @ a @ tm.T


def _rota2(dp: float, ip: float, ar: float) -> np.ndarray:
    """Equivalent de `rota2` (anisotropie.f:930-964) : matrice de rotation a
    partir d'un pole (dp,ip) et d'un angle de rotation propre `ar` (deg)
    autour de ce pole."""
    w = math.pi / 180.0
    tm = np.zeros((3, 3))
    tm[2, :] = cart(1.0, dp, ip)
    tm[0, :] = cart(1.0, dp, ip - 90.0)
    x2 = -tm[0, 1] * tm[2, 2] + tm[2, 1] * tm[0, 2]
    y2 = -tm[2, 0] * tm[0, 2] + tm[2, 2] * tm[0, 0]
    z2 = -tm[0, 0] * tm[2, 1] + tm[2, 0] * tm[0, 1]
    r = math.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
    tm[1, 0], tm[1, 1], tm[1, 2] = -x2 / r, -y2 / r, -z2 / r
    car, sar = math.cos(ar * w), math.sin(ar * w)
    aro = np.array([[car, -sar, 0.0], [sar, car, 0.0], [0.0, 0.0, 1.0]])
    rot = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            rt = 0.0
            for k in range(3):
                for m in range(3):
                    rt += tm[k, i] * tm[m, j] * aro[k, m]
            rot[i, j] = rt
    return rot


def corpen_tensor(a: np.ndarray, dip: float, strike: float) -> np.ndarray:
    """Equivalent de `corpen` pour un TENSEUR (anisotropie.f:888-911) :
    correction de pendage (bedding tilt) via `rota2(strike, 0.0, dip)` puis
    similarite b_ij = sum_km rot_ik rot_jm a_km."""
    rot = _rota2(strike, 0.0, dip)
    return rot @ a @ rot.T


def apply_orientation(m: AMSMeasurement, orientation: int) -> np.ndarray:
    """Renvoie le tenseur de `m` (matrice 3x3) dans le repere demande :
    1=echantillon (brut), 2=in-situ (corfor), 3=corrige du pendage (corfor
    puis corpen) - meme convention 1/2/3 que STARpaleomag_Py (`selection.apply_orientation`)."""
    a = _tensor_matrix(m)
    if orientation == 1:
        return a
    a = corfor_tensor(a, m.cin, m.caz)
    if orientation == 2:
        return a
    if orientation == 3:
        return corpen_tensor(a, m.dip, m.str_)
    raise ValueError(f"invalid orientation: {orientation}")


# Entete du format .pmagani (tabulations) - MEME schema de colonnes que
# STARpaleomag_Py/calcul._PMAGANI_HEADER (les deux applications partagent le
# meme fichier .pmagani, comme elles partageaient deja l'ancien .ANI) -
# demande explicite utilisateur ("can we write the name of the ani
# extension as .pmagani"). PLUS de cin/caz/dip/str ici - demande
# explicite utilisateur ("the codes cin caz dip str are now not needed
# if the pmagani is linked to the prmag file") : ces champs viennent du
# .prmag associe (voir ams_prmag.PrmagSpecimen + la jointure faite dans
# app.ouvrir_prmag_dialog), pas de ce fichier. Un .pmagani ouvert SEUL
# (sans .prmag - voir app.ouvrir_ani_dialog) laisse donc cin/caz/dip/str
# a leurs valeurs par defaut (0.0) : la reorientation (corfor_tensor)
# n'est fiable qu'apres avoir ouvert le .prmag compagnon.
_PMAGANI_HEADER = [
    "specimen", "code2", "etape",
    "k11", "k22", "k33", "k12", "k23", "k13", "s(SI*1e-5)",
    "n_positions", "sigma", "ftest", "ftest12", "ftest23", "quality",
    "field_Am", "export", "info",
]
# TROIS paliers de largeur possibles pour un fichier .pmagani specimen
# (meme principe que _PMAGANI_MEAN_HEADER_LEGACY_LEN/_NO_EXPORT_LEN) :
# le plus ancien format (ni field_Am ni export, info en dernier), un
# format intermediaire (export ajoute, PAS encore field_Am - session
# precedente), et le format actuel (field_Am + export avant info) -
# demande explicite utilisateur ("the field used (in A/m)... The MFK2
# instrument allow measurements at different field values").
_PMAGANI_HEADER_LEGACY_LEN = len(_PMAGANI_HEADER) - 2
_PMAGANI_HEADER_NO_FIELD_LEN = len(_PMAGANI_HEADER) - 1

# `s` en SI*1e-5 (ex. 7261 = 0.07261 SI) - convention Bartington ("since
# the 80s ... measure in 1e-5 SI assuming a volume of 10cc ... could read
# between 0 and 9999"), gardee par l'utilisateur pour comparer Bartington
# et AGICO - demande explicite utilisateur ("We should add in the header
# that s is in SI*1e-5"), meme note que STARpaleomag_Py/calcul._PMAGANI_UNITS_NOTE.
_PMAGANI_UNITS_NOTE = "# s(SI*1e-5): susceptibility in units of 1e-5 SI (Bartington/AGICO convention, e.g. 7261 = 0.07261 SI)\n"

# Les 2 sections du .pmagani - demande explicite utilisateur ("can we
# have the pmagani as pmagres with two blocks one at a specimen level,
# the next one at site level with both blocks having their own set of
# parameters"), MEMES chaines/colonnes exactes que STARpaleomag_Py/calcul.
# _ANI_SPECIMEN_HEADER/_ANI_MEAN_HEADER/_PMAGANI_MEAN_HEADER (les deux
# applications partagent le meme fichier .pmagani). Schema mean = axes
# propres k1>=k2>=k3 (Jelinek) + dec/inc/ellipse de confiance + P/T/L/F/
# Pprim - PAS un duplicata du schema specimen (k11..k13 bruts).
_ANI_SPECIMEN_HEADER = "#specimen tensor results"
_ANI_MEAN_HEADER = "#site mean tensor results"

_PMAGANI_MEAN_HEADER = [
    "site", "code2", "n",
    "k1", "k2", "k3",
    "dec1", "inc1", "dec2", "inc2", "dec3", "inc3",
    "alpha1_1", "alpha2_1", "alpha1_2", "alpha2_2", "alpha1_3", "alpha2_3",
    "P", "T", "L", "F", "Pprim",
    "tilt_correction", "export", "info",
]
# nombre de colonnes d'un fichier ecrit AVANT l'ajout de la colonne
# dediee "tilt_correction" (le plus ancien format - ni tilt_correction
# ni export) - voir _read_pmagani_mean_rows.
_PMAGANI_MEAN_HEADER_LEGACY_LEN = len(_PMAGANI_MEAN_HEADER) - 2
# nombre de colonnes d'un fichier avec tilt_correction mais ecrit AVANT
# l'ajout de "export" (TensorialMeanResult.export) - demande explicite
# utilisateur ("is it possible to export from AMS_py only the data and
# mean tensors that we want to export").
_PMAGANI_MEAN_HEADER_NO_EXPORT_LEN = len(_PMAGANI_MEAN_HEADER) - 1

# iorient (1/2/3, voir amsapp.orientation / radiobuttons "AMS data") <->
# colonne DEDIEE "tilt_correction" de la ligne mean - demande explicite
# utilisateur ("oui ajouter une colonne avant info", suite a "in the
# file pmagani, there is no information about the IS/TC" : l'info
# existait, mais SEULEMENT en texte libre dans `info`, jamais une
# colonne du tableau - une decision deliberee au depart ("encoder ceci
# dans info evite de re-modifier le schema deja valide cote
# STARpaleomag_Py/calcul") revenue ici). Fichiers plus anciens (colonne
# absente, info seul porteur) restent lisibles - voir _mean_row_to_result.
#
# MEME convention que STARpaleomag_Py/calcul._ORIENT_MODE_TAG/_ORIENT_TO_FILE_CODE
# (colonne "IS/TC" de .pmagres) - PAS un enum CE/IS/CP invente separement :
# "Sa"/"IS"/"TC" comme etiquette, et le code stocke est un pourcentage de
# correction de pendage 0-100 au sens MagIC dir_tilt_correction (0=in-situ,
# 100=apres pendage complet ; "1" = coordonnees echantillon, un cas
# degenere sans notion de pourcentage) - demande explicite utilisateur
# ("CP (use TC) and IS are in comment. it should be a variable from 0 to
# 100 (full TC)"). AMS_Py ne produit aujourd'hui que ces 3 valeurs
# discretes (self.orientation n'a que 3 etats, la rotation de pendage
# n'est pas encore fractionnable) mais le CODE STOCKE est bien la
# convention 0-100, pas un enum a 3 lettres - un futur pourcentage
# intermediaire (correction partielle) resterait lisible sans migration.
_ORIENT_MODE_TAG = {1: "Sa", 2: "IS", 3: "TC"}
_ORIENT_TO_FILE_CODE = {1: "1", 2: "0", 3: "100"}
_FILE_CODE_TO_ORIENT = {"1": 1, "0": 2, "100": 3}


def _fmt_pmagani_stat(v: Optional[float]) -> str:
    return "n.d" if v is None else f"{v:.6g}"


def _format_pmagani_specimen_line(m: AMSMeasurement) -> str:
    """Formate UNE ligne de la section specimen (voir _PMAGANI_HEADER) -
    factorise entre les 3 ecrivains ad hoc qui dupliquaient chacun cette
    liste de champs (ams_asc.archive_asc_file, ams_selection.
    import_legacy_ani, et ouvrir_marquer_export_dialog cote app.py qui
    reecrit le fichier entier) - demande explicite utilisateur ("is it
    possible to export from AMS_py only the data and mean tensors that
    we want to export"), pour que la colonne "export" ajoutee ici soit
    ecrite PARTOUT de la meme facon, sans re-dupliquer une 4e fois cette
    liste de 19 champs."""
    info = f'"{m.info}"' if m.info else '""'
    fields = [
        m.id, m.code2, str(m.etape),
        f"{m.k11:.6E}", f"{m.k22:.6E}", f"{m.k33:.6E}",
        f"{m.k12:.6E}", f"{m.k23:.6E}", f"{m.k13:.6E}",
        f"{m.s:.5f}",
        "n.d" if m.n_positions is None else str(m.n_positions),
        _fmt_pmagani_stat(m.sigma), _fmt_pmagani_stat(m.ftest),
        _fmt_pmagani_stat(m.ftest12), _fmt_pmagani_stat(m.ftest23),
        m.quality or "n.d",
        _fmt_pmagani_stat(m.field),
        m.export or "Y",
        info,
    ]
    return "\t".join(fields) + "\n"


def _insert_pmagani_line(path: str, line: str, is_mean: bool) -> None:
    """MEME fonction que STARpaleomag_Py/calcul._insert_pmagani_line (format de
    fichier partage) : une ligne specimen s'insere en fin de section
    specimen (avant l'en-tete mean si elle existe deja) ; une ligne mean
    s'ajoute TOUJOURS en fin de fichier (section creee au besoin)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        lines = [
            "# pmagani v2 - companion of .prmag/.pmagres, join key = specimen "
            "(specimen section) / site (site mean section)",
            _PMAGANI_UNITS_NOTE.rstrip("\n"),
            _ANI_SPECIMEN_HEADER,
            "\t".join(_PMAGANI_HEADER),
        ]
    else:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()

    mean_header_idx = next(
        (i for i, l in enumerate(lines) if l.strip() == _ANI_MEAN_HEADER), None)

    if is_mean:
        if mean_header_idx is None:
            if lines and lines[-1].strip():
                lines.append("")
            lines.append(_ANI_MEAN_HEADER)
            lines.append("\t".join(_PMAGANI_MEAN_HEADER))
        lines.append(line.rstrip("\n"))
    elif mean_header_idx is None:
        lines.append(line.rstrip("\n"))
    else:
        insert_idx = mean_header_idx
        while insert_idx > 0 and not lines[insert_idx - 1].strip():
            insert_idx -= 1
        lines.insert(insert_idx, line.rstrip("\n"))

    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


def _format_pmagani_mean_line(
    site: str, code2: str, result: "TensorialMeanResult", orientation: int, info: str = "",
) -> str:
    """Formate UNE ligne de la section mean pour `result` (3 axes k1>=k2>=
    k3, deja tries - voir ams_stats.TensorialMeanResult). `orientation`
    (1/2/3) ecrite dans la colonne DEDIEE "tilt_correction" (voir
    _ORIENT_TO_FILE_CODE, MEME convention que la colonne "IS/TC" de
    STARpaleomag_Py/calcul .pmagres) - PLUS dupliquee en texte dans
    `info` (l'ancienne convention "tilt_correction: N; ..." restait la
    SEULE source avant l'ajout de cette colonne - demande explicite
    utilisateur "oui ajouter une colonne avant info" ; `info` redevient
    un commentaire libre, plus jamais parse pour une valeur reelle sur
    un fichier ecrit a partir d'ici).
    Leve ValueError si `result` est isotrope/sans axes (rien de pertinent
    a ecrire - voir ams_stats.format_tsmean_box, meme garde)."""
    from ams_stats import shape_params
    if result.ellipsoid_type == 0 or len(result.axes) != 3:
        raise ValueError("isotropic (or incomplete) mean result has no principal axes to save")
    k1, k2, k3 = (ax.eigenvalue for ax in result.axes)
    sp = shape_params(k1, k2, k3)
    tilt_code = _ORIENT_TO_FILE_CODE.get(orientation, "")
    fields = [
        site, code2, str(result.n),
        f"{k1:.6E}", f"{k2:.6E}", f"{k3:.6E}",
        f"{result.axes[0].dec:.2f}", f"{result.axes[0].inc:.2f}",
        f"{result.axes[1].dec:.2f}", f"{result.axes[1].inc:.2f}",
        f"{result.axes[2].dec:.2f}", f"{result.axes[2].inc:.2f}",
        f"{result.axes[0].alpha[0]:.3f}", f"{result.axes[0].alpha[1]:.3f}",
        f"{result.axes[1].alpha[0]:.3f}", f"{result.axes[1].alpha[1]:.3f}",
        f"{result.axes[2].alpha[0]:.3f}", f"{result.axes[2].alpha[1]:.3f}",
        _fmt_pmagani_stat(sp["P"]), _fmt_pmagani_stat(sp["T"]),
        _fmt_pmagani_stat(sp["L"]), _fmt_pmagani_stat(sp["F"]), _fmt_pmagani_stat(sp["Pprim"]),
        tilt_code,
        result.export or "Y",
        f'"{info}"' if info else '""',
    ]
    return "\t".join(fields) + "\n"


def write_ani_mean_result(
    path: str, site: str, code2: str, result: "TensorialMeanResult",
    orientation: int, info: str = "",
) -> None:
    """Ajoute UNE moyenne de site a la section '#site mean tensor
    results' de `path` (creee au besoin - voir _insert_pmagani_line) -
    demande explicite utilisateur ("on import, the results are not saved
    in the file pmagani"). Leve ValueError si `result` est isotrope (voir
    _format_pmagani_mean_line)."""
    _insert_pmagani_line(path, _format_pmagani_mean_line(site, code2, result, orientation, info), is_mean=True)


def _read_pmagani_mean_rows(path: str) -> List[dict]:
    """Lit la section '#site mean tensor results' de `path` (voir
    _PMAGANI_MEAN_HEADER) en lignes brutes (dict) - [] si absente."""
    rows: List[dict] = []
    in_mean = False
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if stripped == _ANI_MEAN_HEADER:
                    in_mean = True
                continue
            if not in_mean:
                continue
            parts = line.split("\t")
            if not parts or parts[0] == "site":
                continue  # ligne d'entete
            # Retro-compatibilite sur le nombre de colonnes (meme
            # principe que calcul._read_pmagani_mean_tensors, meme
            # fichier partage) - 3 paliers : le plus ancien format (ni
            # tilt_correction ni export, info en dernier), le format
            # intermediaire (tilt_correction ajoute, PAS encore export),
            # et le format actuel (tilt_correction + export avant info).
            if len(parts) >= len(_PMAGANI_MEAN_HEADER):
                tilt_correction, export_raw, info_raw = parts[23].strip(), parts[24].strip(), parts[25]
            elif len(parts) >= _PMAGANI_MEAN_HEADER_NO_EXPORT_LEN:
                tilt_correction, export_raw, info_raw = parts[23].strip(), "Y", parts[24]
            elif len(parts) >= _PMAGANI_MEAN_HEADER_LEGACY_LEN:
                tilt_correction, export_raw, info_raw = "", "Y", parts[23]
            else:
                continue
            try:
                rows.append({
                    "site": parts[0], "code2": parts[1], "n": int(float(parts[2])),
                    "k1": float(parts[3]), "k2": float(parts[4]), "k3": float(parts[5]),
                    "dec1": float(parts[6]), "inc1": float(parts[7]),
                    "dec2": float(parts[8]), "inc2": float(parts[9]),
                    "dec3": float(parts[10]), "inc3": float(parts[11]),
                    "alpha1_1": float(parts[12]), "alpha2_1": float(parts[13]),
                    "alpha1_2": float(parts[14]), "alpha2_2": float(parts[15]),
                    "alpha1_3": float(parts[16]), "alpha2_3": float(parts[17]),
                    "tilt_correction": tilt_correction,
                    "export": export_raw if export_raw in ("Y", "N") else "Y",
                    "info": info_raw.strip('"'),
                })
            except ValueError:
                continue
    return rows


def _restm_from_alpha(alpha1_deg: float, alpha2_deg: float) -> "np.ndarray":
    """Reconstruit une `restm` (2x2, voir ams_stats.ellips/AxisResult)
    depuis les seuls demi-angles alpha1(grand axe)/alpha2(petit axe) en
    degres DEJA sauvegardes (.pmagani ne stocke pas `v`, l'orientation
    propre 2x2 de l'ellipse dans son plan tangent) - demande explicite
    utilisateur ("is it possible to... plot the ellipses") : sans ceci,
    une moyenne RELUE depuis un .pmagani (par opposition a une moyenne
    FRAICHEMENT calculee dans la meme session, ou `ellips()` peuple
    `restm` directement) n'affiche jamais d'ellipse - draw_stereo_
    confidence_ellipses saute silencieusement tout axe a `restm is
    None`. Assume l'ellipse ALIGNEE sur les deux autres axes propres
    (grand axe vers l'axe j, petit axe vers l'axe k) - PAS une
    approximation nouvelle, c'est deja la simplification documentee et
    appliquee ailleurs dans ce module pour eta/zeta (voir ams_stats.
    magic_site_aniso_fields : "eta/zeta pointent vers les DEUX AUTRES
    axes propres eux-memes, pas une orientation d'ellipse recalculee").
    `ellips()` donne alpha=atan(sqrt(c*lambda)) et restm=v@diag(sqrt(c*
    lambda))@v.T ; v=Identite ici donne restm=diag(tan(alpha1),
    tan(alpha2)) - la vraie rotation `v` (si l'ellipse n'etait pas deja
    axis-aligned dans le fichier d'origine) reste perdue, l'ellipse
    redessinee est donc une approximation axis-aligned, pas un
    round-trip exact."""
    return np.diag([
        math.tan(math.radians(alpha1_deg)),
        math.tan(math.radians(alpha2_deg)),
    ])


def _mean_row_to_result(row: dict) -> Tuple[int, "TensorialMeanResult"]:
    """Convertit une ligne brute (voir _read_pmagani_mean_rows) en
    (orientation, TensorialMeanResult) - meme convention de tuple que
    app.py:self.mean_results. `orientation` retrouvee en priorite dans
    la colonne dediee `tilt_correction` ; a defaut (fichier ecrit avant
    son ajout), repli sur le texte libre `info` (ancienne convention
    "tilt_correction: N; ..."). 2 (in situ, code "0") si ni l'une ni
    l'autre n'est exploitable."""
    import re
    from ams_stats import AxisResult, TensorialMeanResult
    orientation = 2
    tilt_correction = (row.get("tilt_correction") or "").strip()
    if tilt_correction in _FILE_CODE_TO_ORIENT:
        orientation = _FILE_CODE_TO_ORIENT[tilt_correction]
    else:
        m = re.search(r"tilt_correction:\s*(1|0|100)\b", row.get("info", ""))
        if m:
            orientation = _FILE_CODE_TO_ORIENT.get(m.group(1), 2)
    axes = [
        AxisResult(eigenvalue=row["k1"], dec=row["dec1"], inc=row["inc1"],
                   alpha=(row["alpha1_1"], row["alpha2_1"]),
                   restm=_restm_from_alpha(row["alpha1_1"], row["alpha2_1"])),
        AxisResult(eigenvalue=row["k2"], dec=row["dec2"], inc=row["inc2"],
                   alpha=(row["alpha1_2"], row["alpha2_2"]),
                   restm=_restm_from_alpha(row["alpha1_2"], row["alpha2_2"])),
        AxisResult(eigenvalue=row["k3"], dec=row["dec3"], inc=row["inc3"],
                   alpha=(row["alpha1_3"], row["alpha2_3"]),
                   restm=_restm_from_alpha(row["alpha1_3"], row["alpha2_3"])),
    ]
    result = TensorialMeanResult(
        n=row["n"], axes=axes, ellipsoid_type=1, id=row["site"],
        export=row.get("export", "Y"))
    return orientation, result


def read_ani_mean_results_from_pmagani(path: str) -> List[Tuple[int, "TensorialMeanResult", str]]:
    """Lit la section mean d'un .pmagani deja ecrite (par write_ani_mean_result
    - cote AMS_Py ou STARpaleomag_Py, meme fichier partage) et la reconvertit en
    (orientation, TensorialMeanResult, code2) - le 3eme element (code2,
    absent de app.py:self.mean_results) est retourne a part pour laisser
    l'appelant filtrer/afficher sans perte d'info."""
    out = []
    for row in _read_pmagani_mean_rows(path):
        orientation, result = _mean_row_to_result(row)
        out.append((orientation, result, row["code2"]))
    return out


def mark_pmagani_export(
    path: str, selected_specimen_ids: "set", selected_mean_keys: "set",
) -> Tuple[int, int, int, int]:
    """Reecrit .pmagani EN PLACE, colonne "export" UNIQUEMENT (Y pour les
    lignes selectionnees, N pour toutes les AUTRES lignes du MEME
    fichier) - demande explicite utilisateur ("is it possible to export
    from AMS_py only the data and mean tensors that we want to export
    and only these selected data will be taken into account in the main
    export from Starpaleomag"). TOUTES les autres colonnes DEJA presentes
    sur chaque ligne restent BYTE POUR BYTE identiques (edition directe
    de la colonne "export" dans la liste `parts` deja tabulee, PAS un
    reformatage complet depuis des objets reconstruits - evite tout
    risque de deriver la precision/le format d'une colonne non
    concernee, ex. P/T/L/F/Pprim qui ne sont meme pas conservees par
    _read_pmagani_mean_rows).

    `selected_specimen_ids` : ensemble d'id de specimen (MAJUSCULES) a
    marquer export=Y - tout specimen DU FICHIER absent de cet ensemble
    est marque export=N. `selected_mean_keys` : ensemble de
    (site, code2, tilt_correction_code_fichier) en MAJUSCULES pour
    site/code2 (voir _ORIENT_TO_FILE_CODE pour convertir une orientation
    1/2/3 en code fichier "1"/"0"/"100"), meme principe pour le reste.

    Une ligne "format ancien" (ecrite avant l'ajout de la colonne
    "export" - et pour la section mean, potentiellement avant
    "tilt_correction" aussi) n'a PAS ces colonnes du tout : elles sont
    INSEREES (pas seulement editees) a la bonne position - premiere
    reecriture reelle d'un vrai fichier utilisateur qui predate cette
    fonctionnalite, sinon la colonne "export" ne serait JAMAIS visible
    dans un tel fichier (le point mort observe : "export" absent meme
    apres avoir lance "Mark selection for MagIC export..."). Une ligne
    mean sans "tilt_correction" du tout n'a pas non plus de code fichier
    connu pour la cle de correspondance - `tilt_correction_code_fichier`
    vaut alors "" (ne correspondra qu'a une cle appelante elle-meme
    construite avec un code "", ce qui n'arrive normalement jamais - une
    telle ligne tres ancienne est donc marquee N faute de mieux). Les
    lignes d'en-tete ("specimen ..."/"site ...") sont elles aussi mises a
    niveau vers _PMAGANI_HEADER/_PMAGANI_MEAN_HEADER si necessaire (noms
    de colonnes uniquement, aucune donnee). Une ligne encore plus courte
    (format tres ancien, en dessous meme du palier legacy le plus bas) est
    laissee TELLE QUELLE (ni touchee ni comptee) - reste lue "export=Y"
    par defaut par les deux applications.

    Retourne (specimens_inclus, specimens_exclus, moyennes_incluses,
    moyennes_exclues)."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    in_mean = False
    n_spec_in = n_spec_out = n_mean_in = n_mean_out = 0
    out_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("#"):
            if stripped == _ANI_MEAN_HEADER:
                in_mean = True
            out_lines.append(line)
            continue
        if not stripped:
            out_lines.append(line)
            continue
        parts = line.split("\t")
        if not in_mean:
            if not parts:
                out_lines.append(line)
                continue
            if parts[0] == "specimen":
                out_lines.append("\t".join(_PMAGANI_HEADER) if len(parts) < len(_PMAGANI_HEADER) else line)
                continue
            if len(parts) >= len(_PMAGANI_HEADER):
                pass  # deja field_Am + export en place, index 16/17
            elif len(parts) >= _PMAGANI_HEADER_NO_FIELD_LEN:
                parts = parts[:16] + ["n.d"] + parts[16:]  # insere "field_Am" avant "export"
            elif len(parts) >= _PMAGANI_HEADER_LEGACY_LEN:
                parts = parts[:16] + ["n.d", ""] + parts[16:]  # insere field_Am="n.d" puis export avant "info"
            else:
                out_lines.append(line)
                continue
            included = parts[0].strip().upper() in selected_specimen_ids
            parts[17] = "Y" if included else "N"
            if included:
                n_spec_in += 1
            else:
                n_spec_out += 1
            out_lines.append("\t".join(parts))
        else:
            if not parts:
                out_lines.append(line)
                continue
            if parts[0] == "site":
                out_lines.append("\t".join(_PMAGANI_MEAN_HEADER) if len(parts) < len(_PMAGANI_MEAN_HEADER) else line)
                continue
            if len(parts) >= len(_PMAGANI_MEAN_HEADER):
                pass  # deja tilt_correction+export en place, index 23/24
            elif len(parts) >= _PMAGANI_MEAN_HEADER_NO_EXPORT_LEN:
                parts = parts[:24] + [""] + parts[24:]  # insere "export" avant "info"
            elif len(parts) >= _PMAGANI_MEAN_HEADER_LEGACY_LEN:
                parts = parts[:23] + ["", ""] + parts[23:]  # insere tilt_correction="" puis "export" avant "info"
            else:
                out_lines.append(line)
                continue
            key = (parts[0].strip().upper(), parts[1].strip().upper(), parts[23].strip())
            included = key in selected_mean_keys
            parts[24] = "Y" if included else "N"
            if included:
                n_mean_in += 1
            else:
                n_mean_out += 1
            out_lines.append("\t".join(parts))

    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(out_lines) + "\n")
    return n_spec_in, n_spec_out, n_mean_in, n_mean_out


def _parse_pmagani_stat(v: str) -> Optional[float]:
    v = v.strip()
    if not v or v == "n.d":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _read_pmagani_file(path: str) -> List[AMSMeasurement]:
    """Lit un fichier .pmagani (nouveau format tabule avec entete, voir
    STARpaleomag_Py/calcul._read_pmagani_tensors - meme format de colonnes).
    cin/caz/dip/str_ restent a leur valeur par defaut (0.0) - PLUS dans ce
    fichier (voir commentaire au-dessus de _PMAGANI_HEADER) ; a completer
    depuis le .prmag associe si besoin (voir app.ouvrir_prmag_dialog).

    S'ARRETE des que l'en-tete '#site mean tensor results' est rencontre
    (voir _ANI_MEAN_HEADER) : cette 2e section (moyennes de site, schema
    Jelinek different) n'est pas encore geree cote AMS_Py et ne doit
    jamais etre lue comme des tenseurs specimen."""
    out: List[AMSMeasurement] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if stripped == _ANI_MEAN_HEADER:
                    break
                continue
            parts = line.split("\t")
            if not parts or parts[0] == "specimen":
                continue  # ligne d'entete
            if len(parts) < 9:
                continue
            try:
                m = AMSMeasurement(
                    id=parts[0], code2=parts[1].ljust(2)[:2],
                    etape=int(float(parts[2])),
                    k11=float(parts[3]), k22=float(parts[4]), k33=float(parts[5]),
                    k12=float(parts[6]), k23=float(parts[7]), k13=float(parts[8]),
                    s=float(parts[9]) if len(parts) > 9 else 0.0,
                )
            except (ValueError, IndexError):
                continue
            if len(parts) > 10:
                try:
                    m.n_positions = int(float(parts[10]))
                except ValueError:
                    m.n_positions = None
            m.sigma = _parse_pmagani_stat(parts[11]) if len(parts) > 11 else None
            m.ftest = _parse_pmagani_stat(parts[12]) if len(parts) > 12 else None
            m.ftest12 = _parse_pmagani_stat(parts[13]) if len(parts) > 13 else None
            m.ftest23 = _parse_pmagani_stat(parts[14]) if len(parts) > 14 else None
            if len(parts) > 15 and parts[15].strip() in ("g", "b"):
                m.quality = parts[15].strip()
            # field_Am/export/info : retro-compatibilite sur le nombre de
            # colonnes (meme principe que _PMAGANI_MEAN_HEADER_LEGACY_LEN/
            # _NO_EXPORT_LEN) - TROIS paliers : le plus ancien fichier n'a
            # ni field_Am ni export (info en 16), un fichier intermediaire
            # a export mais pas field_Am (export en 16, info en 17), le
            # format actuel a field_Am puis export avant info (16/17/18).
            if len(parts) >= len(_PMAGANI_HEADER):
                field_raw, export_raw, info_idx = parts[16].strip(), parts[17].strip(), 18
            elif len(parts) >= _PMAGANI_HEADER_NO_FIELD_LEN:
                field_raw, export_raw, info_idx = "", parts[16].strip(), 17
            elif len(parts) >= _PMAGANI_HEADER_LEGACY_LEN:
                field_raw, export_raw, info_idx = "", "Y", 16
            else:
                field_raw, export_raw, info_idx = "", "Y", None
            m.field = _parse_pmagani_stat(field_raw)
            m.export = export_raw if export_raw in ("Y", "N") else "Y"
            if info_idx is not None and len(parts) > info_idx:
                info = parts[info_idx].strip()
                if info.startswith('"') and info.endswith('"') and len(info) >= 2:
                    info = info[1:-1]
                m.info = info
            out.append(m)
    return out


def _read_ani_file_legacy(path: str) -> List[AMSMeasurement]:
    """Equivalent de `openani` (lect_asc.f:592-676) : lit un ANCIEN fichier
    .ANI AMS (lignes 'D ...' = mesures brutes par position ; les lignes
    'R ...' = resultats tensoriels deja calcules sont ignorees ICI, voir
    `read_ani_results` - fonction separee, meme fichier). Format D
    verifie octet-pres contre un fichier .ANI reel
    (reference/AMS_OSX_AWE/NANGQIAN.ASC.ANI) : `D id cin caz dip str
    etape code2 k11 k22 k33 k12 k23 k13 s ["info"]` - le champ info
    est optionnel (branche de secours 1101 du Fortran en son absence,
    comme observe sur ce fichier reel)."""
    out: List[AMSMeasurement] = []
    with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
        for line in f:
            parts = line.split(None, 14)
            if not parts or parts[0] != "D":
                continue
            if len(parts) < 14:
                continue
            try:
                m = AMSMeasurement(
                    id=parts[1],
                    cin=float(parts[2]), caz=float(parts[3]),
                    dip=float(parts[4]), str_=float(parts[5]),
                    etape=int(float(parts[6])), code2=parts[7].ljust(2)[:2],
                    k11=float(parts[8]), k22=float(parts[9]), k33=float(parts[10]),
                    k12=float(parts[11]), k23=float(parts[12]), k13=float(parts[13]),
                    s=float(parts[14].split()[0]) if len(parts) > 14 else 0.0,
                )
            except (ValueError, IndexError):
                continue
            if len(parts) > 14:
                rest = parts[14]
                q1 = rest.find('"')
                if q1 >= 0:
                    q2 = rest.find('"', q1 + 1)
                    m.info = rest[q1 + 1:q2] if q2 > q1 else ""
            out.append(m)
    return out


def mean_result_site(id_: str) -> Optional[str]:
    """Extrait le nom de site d'un `id` de resultat de moyenne tensorielle,
    quand c'est possible - port de la convention de `storeres`
    (anisotropie.f:3041, JAMAIS encore ecrite par AMS_Py lui-meme - voir
    read_ani_mean_results) : `id` concatene les 6 premiers caracteres du
    PREMIER et du DERNIER specimen de la liste editee au moment du calcul
    (ex. "12TU7412TU74"). Quand les deux moities de 6 caracteres sont
    identiques, tous les specimens moyennes partageaient le meme prefixe
    (convention 6 caracteres = annee(2)+site(4), meme convention que
    STARpaleomag_Py/field_notes.py) - demande explicite utilisateur
    ("12TU7412TU74 indicate that only samples from site 12TU74 were used
    for the mean") : retourne ce prefixe. Sinon (moities differentes, ou
    `id` pas exactement 12 caracteres), la moyenne melange plusieurs
    sites/specimens et aucun site unique n'en est extractible - retourne
    None plutot que de deviner."""
    if len(id_) != 12:
        return None
    first, last = id_[:6], id_[6:]
    return first if first == last else None


def read_ani_mean_results(path: str) -> List[Tuple[int, "TensorialMeanResult"]]:
    """Lit les lignes 'R ...' (resultats de moyenne tensorielle DEJA
    calcules, ex. par une ancienne version AWE de ce logiciel) d'un
    fichier de resultats ANI legacy - port de la lecture de `storeres`/
    `sauveresfich` (anisotropie.f:2989-3049/3210-3241) - demande explicite
    utilisateur ("Is it possible to import the line with a mean tensor
    (start with R)"). PAS encore ecrites par AMS_Py lui-meme (le
    'save results' equivalent - sauveresfich - n'est pas encore porte,
    voir ouvrir_tsmean_dialog) : cette fonction permet neanmoins
    d'IMPORTER des resultats deja calcules et sauvegardes par l'ancienne
    application Fortran.

    Format (list-directed Fortran, `write(27,*) "R ", id, n, iorient,
    ...` - un tenseur moyen `TensorialMeanResult` complet a 3 axes) :
        R <id> <n> <iorient>
          <k1> <#erreur k1> <dec k1> <inc k1>
          <k2> <#erreur k2> <dec k2> <inc k2>
          <k3> <#erreur k3> <dec k3> <inc k3>
          <tmres(1,1,1)> <tmres(2,1,1)> <tmres(3,1,1)>
          <tmres(1,2,1)> <tmres(2,2,1)> <tmres(3,2,1)>
          <tmres(1,1,2)> <tmres(2,1,2)> <tmres(3,1,2)>
          <tmres(1,2,2)> <tmres(2,2,2)> <tmres(3,2,2)>
          <alpha1 k1> <alpha1 k2> <alpha1 k3>
          <alpha2 k1> <alpha2 k2> <alpha2 k3>
          <kmax> <kint> <kmin>
    (37 tokens au total sur la ligne : "R", id, n, iorient, puis 33
    valeurs numeriques = 12 (axk1) + 12 (tmres) + 6 (alphares) + 3
    (kmax/kint/kmin) - voir le decoupage exact ci-dessous). Ordre de
    tmres/alphares VERIFIE par symetrie sur un exemple reel fourni par
    l'utilisateur (12TU74, tmres[j,2,1]==tmres[j,1,2] pour les 3 axes,
    confirmant l'entrelacement j-plus-rapide/k/l issu du triple do-loop
    Fortran) et par kmax/kint/kmin == les 3 valeurs propres k1/k2/k3.
    `kmax`/`kint`/`kmin` (redondants avec les valeurs propres, gardes par
    le Fortran pour un acces direct) ne sont PAS reportes separement dans
    `TensorialMeanResult` (deja dans `axes[i].eigenvalue`) - juste
    verifies pour rejeter une ligne visiblement corrompue.

    `mean_tensor_normalized`/`ellipsoid_type` de `TensorialMeanResult` ne
    sont PAS reconstructibles depuis ce format (le tenseur moyen brut
    n'y est pas ecrit, seule sa decomposition propre) - laisses a leurs
    valeurs par defaut (None/-1). Le site d'un resultat mono-site peut
    ensuite etre retrouve via `mean_result_site(result.id)`.

    Retourne une liste de (iorient, TensorialMeanResult) - meme
    convention de tuple que app.py:self.mean_results, pour un ajout
    direct sans transformation supplementaire."""
    from ams_stats import AxisResult, TensorialMeanResult

    out: List[Tuple[int, TensorialMeanResult]] = []
    with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
        for line in f:
            parts = line.split()
            if not parts or parts[0] != "R":
                continue
            if len(parts) < 37:
                continue
            try:
                id_ = parts[1]
                n = int(float(parts[2]))
                iorient = int(float(parts[3]))
                nums = [float(x) for x in parts[4:37]]
            except (ValueError, IndexError):
                continue

            axk1 = nums[0:12]
            tmres_flat = nums[12:24]
            alphares_flat = nums[24:30]
            kmax, kint, kmin = nums[30:33]

            axes: List[AxisResult] = []
            for j in range(3):
                eigenvalue, std, dec, inc = axk1[j * 4:j * 4 + 4]
                # voir docstring : tmres_flat[0+j]=tmres(j+1,1,1),
                # [3+j]=tmres(j+1,2,1), [6+j]=tmres(j+1,1,2), [9+j]=tmres(j+1,2,2)
                restm = np.array([
                    [tmres_flat[0 + j], tmres_flat[6 + j]],
                    [tmres_flat[3 + j], tmres_flat[9 + j]],
                ])
                alpha = (alphares_flat[j], alphares_flat[3 + j])
                axes.append(AxisResult(
                    eigenvalue=eigenvalue, std=std, dec=dec, inc=inc,
                    restm=restm, alpha=alpha,
                ))

            # kmax/kint/kmin redondants avec axes[0/1/2].eigenvalue - juste
            # une verification de coherence, pas une donnee supplementaire.
            if abs(kmax - axes[0].eigenvalue) > 1e-3 or abs(kmin - axes[2].eigenvalue) > 1e-3:
                continue

            out.append((iorient, TensorialMeanResult(n=n, axes=axes, id=id_)))
    return out


def read_ani_file(path: str) -> List[AMSMeasurement]:
    """Lit un fichier .pmagani (nouveau format tabule, voir
    _read_pmagani_file) OU un ancien .ANI (format D, voir
    _read_ani_file_legacy) - dispatch par EXTENSION : les deux formats
    restent lisibles - demande explicite utilisateur ("can we import old
    style .ANI in these pmagani style"). Pour convertir un ancien fichier
    vers le nouveau format, voir import_legacy_ani."""
    if not os.path.exists(path):
        return []
    ext = os.path.splitext(path)[1].lower()
    if ext == ".pmagani":
        return _read_pmagani_file(path)
    return _read_ani_file_legacy(path)


def import_legacy_ani(old_path: str, new_path: Optional[str] = None) -> str:
    """Convertit un ancien fichier .ANI (D id cin caz dip str etape code2
    k11 k22 k33 k12 k23 k13 s [info]) vers le nouveau format .pmagani
    (tabulations, entete de colonnes) - demande explicite utilisateur
    ("can we import old style .ANI in these pmagani style"), equivalent
    de STARpaleomag_Py/calcul.import_legacy_ani (meme format de sortie). Les
    colonnes nouvelles (n_positions/sigma/ftest/ftest12/ftest23) sont
    laissees "n.d" : non disponibles dans l'ancien format. `new_path` par
    defaut : voir ams_prmag.ani_path_for(old_path).

    Reporte aussi les lignes 'R ...' (moyennes de site DEJA calculees par
    l'ancienne application, voir read_ani_mean_results) dans la section
    '#site mean tensor results' du nouveau .pmagani - demande explicite
    utilisateur ("can we check the new import legacy files where the IS
    TC will also be written from 0 to 100") : read_ani_mean_results
    existait deja (demande anterieure "Is it possible to import the line
    with a mean tensor") mais n'etait jusqu'ici appelee nulle part - ces
    resultats etaient donc SILENCIEUSEMENT PERDUS a l'import, jamais
    ecrits dans le nouveau fichier. Chaque resultat porte deja son propre
    `iorient` d'origine (Sample/IS/TC), encode ici avec la MEME convention
    0-100 que les moyennes fraichement calculees (voir
    _format_pmagani_mean_line/_ORIENT_TO_FILE_CODE) - pas de code2 dans le
    format R d'origine, "N0" (AMS) par defaut. Les resultats iorient==4
    ("Dr", sans equivalent dans le modele a 3 etats d'AMS_Py -
    self.orientation) sont ECARTES plutot que mal-etiquetes IS par
    defaut (voir _mean_row_to_result, dont c'est precisement le repli
    silencieux qui motive de ne PAS y laisser tomber un code non
    reconnu)."""
    from ams_prmag import ani_path_for
    if new_path is None:
        new_path = ani_path_for(old_path)
    measurements = _read_ani_file_legacy(old_path)
    with open(new_path, "w", encoding="utf-8") as out:
        out.write(f"# pmagani v1 - imported from legacy .ANI: {os.path.basename(old_path)}\n")
        out.write(_PMAGANI_UNITS_NOTE)
        out.write("\t".join(_PMAGANI_HEADER) + "\n")
        for m in measurements:
            out.write(_format_pmagani_specimen_line(m))

    for iorient, result in read_ani_mean_results(old_path):
        if iorient not in _ORIENT_TO_FILE_CODE:
            continue  # iorient==4 ("Dr") : pas de code 0-100 fiable, voir docstring
        site = mean_result_site(result.id) or result.id.strip()
        try:
            write_ani_mean_result(
                new_path, site, "N0", result, iorient,
                info=f"imported from legacy .ANI: {os.path.basename(old_path)}")
        except ValueError:
            continue  # tenseur isotrope/incomplet - rien de pertinent a ecrire
    return new_path


def select_measurements(
    measurements: List[AMSMeasurement],
    pattern: str = "*",
    step_min: int = 0,
    step_max: int = 9000,
    demag1: str = "*",
    demag2: str = "*",
) -> List[AMSMeasurement]:
    """Equivalent de `selmes` (lect_asc.f:436-561) : filtre `measurements`
    par prefixe d'identifiant (nlen premiers caracteres de `pattern`, '*' =
    tout), bornes de `etape`, et code2 (1er/2e caractere, '*' = tout) -
    meme convention que `selection.select_samples` cote STARpaleomag_Py."""
    pattern = (pattern or "*").upper().strip()[:12] or "*"
    nlen = len(pattern)
    out = []
    for m in measurements:
        mid = (m.id or "").ljust(12).upper()
        if pattern[0] != "*" and mid[:nlen] != pattern[:nlen]:
            continue
        if not (step_min <= m.etape <= step_max):
            continue
        c2 = (m.code2 or "  ").ljust(2)
        if demag1 != "*" and c2[0:1] != demag1:
            continue
        if demag2 != "*" and c2[1:2] != demag2:
            continue
        out.append(m)
    return out


def list_measurements(measurements: List[AMSMeasurement], out=None) -> None:
    """Equivalent de `liste`/`lismes` (anisotropie.f:1250, lect_asc.f:220-226) :
    affichage tabulaire des mesures brutes."""
    import sys
    import io as _io
    if out is None:
        out = sys.stdout
    for i, m in enumerate(measurements, start=1):
        out.write(
            f"{i:4d}  {m.id:<12} {m.cin:6.1f} {m.caz:6.1f} {m.dip:6.1f} {m.str_:6.1f} "
            f"{m.etape:5d} {m.code2:2s}  "
            f"{m.k11:9.5f} {m.k22:9.5f} {m.k33:9.5f} {m.k12:10.3E} {m.k23:10.3E} {m.k13:10.3E}  "
            f"s={m.s:9.4f}"
            + (f"  {m.info}" if m.info else "")
            + "\n"
        )
