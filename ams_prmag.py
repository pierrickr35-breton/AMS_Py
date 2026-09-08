"""Lecture minimale des fichiers .prmag (format STARpaleomag_Py) pour AMS_Py.

AMS_Py ne traite jamais les mesures elles-memes (celles-ci vivent dans le
.ANI, deja au format tenseur - voir ams_selection.read_ani_file) : le seul
besoin est de recuperer, par specimen, le site/sample MagIC declares dans le
fichier .prmag associe, pour enrichir l'export MagIC ("in that case, it will
be much easier to export Anisotropy data to Magic" - demande explicite
utilisateur). Seule la premiere ligne d'entete de chaque bloc (specimen/
sample/site/...) est donc lue ; le reste du bloc (orientation, stratigraphie,
mesures) est ignore.

Le fichier .ANI compagnon garde volontairement seulement le numero de
specimen ("we just keep the specimen number in .ani" - demande explicite
utilisateur) : c'est ce module qui fait le lien vers le site/sample complet,
par jointure sur ce numero de specimen - meme convention de nommage que
STARpaleomag_Py/calcul.ani_path_for (memes basenames pour .prmag/.pmagres/.ANI)."""

import os
from dataclasses import dataclass
from typing import Dict, Optional


def ani_path_for(data_path: str) -> str:
    """Equivalent de calcul.ani_path_for (STARpaleomag_Py) : le fichier .pmagani
    compagnon d'un .prmag/.pmagres partage le meme nom de base. Extension
    ".pmagani" (et non plus ".ANI") - demande explicite utilisateur ("can
    we write the name of the ani extension as .pmagani") ; les anciens
    .ANI restent lisibles directement (voir read_ani_file dans
    ams_selection.py, dispatch par extension) ou peuvent etre convertis
    (voir ams_selection.import_legacy_ani)."""
    base, _ext = os.path.splitext(data_path)
    return base + ".pmagani"


@dataclass
class PrmagSpecimen:
    id: str = ""
    site: str = ""
    sample: str = ""
    # cin/caz (orientation du carottier) et dip/str_ (pendage/direction de
    # stratification) - MEMES champs, memes noms, que testlect.Pmag (voir
    # read_prmag_file cote STARpaleomag_Py) : .pmagani ne les stocke plus lui-
    # meme ("the codes cin caz dip str are now not needed if the pmagani
    # is linked to the prmag file" - demande explicite utilisateur), ce
    # sont maintenant CES valeurs-la (jointure par specimen, voir
    # app.ouvrir_prmag_dialog) qui alimentent AMSMeasurement.cin/caz/dip/
    # str_ pour la reorientation (corfor_tensor).
    cin: float = 0.0
    caz: float = 0.0
    dip: float = 0.0
    str_: float = 0.0


def _prmag_kv_line(line: str) -> dict:
    result = {}
    for chunk in line.split("\t"):
        if ":" not in chunk:
            continue
        k, _sep, v = chunk.partition(":")
        result[k.strip()] = v.strip()
    return result


def _prmag_text(v: Optional[str]) -> str:
    v = (v or "").strip()
    return "" if v == "n.d" else v


def _prmag_nd(txt: Optional[str], default: float = 0.0) -> float:
    txt = (txt or "").strip()
    if not txt or txt == "n.d":
        return default
    try:
        return float(txt)
    except ValueError:
        return default


def read_prmag_specimens(filepath: str, encoding: str = "utf-8") -> Dict[str, PrmagSpecimen]:
    """Lit un fichier .prmag et retourne {specimen_id: PrmagSpecimen}, en ne
    parsant que les 3 premieres lignes d'entete de chaque bloc (specimen/
    sample/site, puis azimuth/dip, puis bed_dip_strike/bed_dip - voir
    testlect.read_prmag_file cote STARpaleomag_Py pour le detail complet du
    format) ; le reste du bloc (4e ligne d'entete, mesures) est saute
    jusqu'a la ligne blanche suivante, sans etre interprete."""
    with open(filepath, "r", encoding=encoding) as f:
        lines = [raw.rstrip("\n") for raw in f]

    specimens: Dict[str, PrmagSpecimen] = {}
    i, n = 0, len(lines)

    def _skip_blank_and_comments():
        nonlocal i
        while i < n and (lines[i].strip() == "" or lines[i].lstrip().startswith("#")):
            i += 1

    _skip_blank_and_comments()
    while i < n:
        line_a = _prmag_kv_line(lines[i])
        line_b = _prmag_kv_line(lines[i + 1]) if i + 1 < n else {}
        line_c = _prmag_kv_line(lines[i + 2]) if i + 2 < n else {}
        specimen_id = _prmag_text(line_a.get("specimen"))
        if specimen_id:
            specimens[specimen_id] = PrmagSpecimen(
                id=specimen_id,
                sample=_prmag_text(line_a.get("sample")),
                site=_prmag_text(line_a.get("site")),
                caz=_prmag_nd(line_b.get("azimuth")),
                cin=_prmag_nd(line_b.get("dip")),
                str_=_prmag_nd(line_c.get("bed_dip_strike")),
                dip=_prmag_nd(line_c.get("bed_dip")),
            )
        while i < n and lines[i].strip() != "":
            i += 1
        _skip_blank_and_comments()

    return specimens
