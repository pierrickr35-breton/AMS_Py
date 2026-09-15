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
import re
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


def _legacy_ani_specimen_to_sample(specimen: str) -> str:
    """19DN1511A -> 19DN1511 (retire la lettre de sous-carotte finale,
    meme convention que detailed_export._sample_display_name cote
    STARpaleomag_Py)."""
    return re.sub(r"[A-Za-z]$", "", specimen.strip())


def _legacy_ani_sample_to_site(sample: str) -> str:
    """19DN1511 -> 19DN15 (les 6 premiers caracteres = annee+site,
    convention Rennes)."""
    return sample[:6]


def create_prmag_from_legacy_ani(
    ani_path: str,
    prmag_path: str,
    volume: float = 10.80,
) -> int:
    """Cree un .prmag STARpaleomag_Py (mesures VIDES) a partir des
    informations d'orientation d'un ANCIEN fichier .ANI - demande
    explicite utilisateur ("un collegue souhaite avoir la possibilite
    de creer le prmag a partir du .ANI qui contient les infos de
    corrections de carotte"). Modifications initiales fournies par un
    collegue (voir historique de conversation), integrees ici avec un
    correctif de deduplication (voir plus bas).

    Colonnes du .ANI (verifie octet-pres contre de vrais fichiers .ANI
    reels, ex. Chili_Briques_pmag.ANI - meme format documente par
    ams_selection._read_ani_file_legacy : `D id cin caz dip str etape
    code2 k11 k22 k33 k12 k23 k13 s [info]`) :
        1 (0-index) = specimen
        2 = dip (cin - pendage/plunge de l'axe du carottier)
        3 = azimuth (caz - azimut de l'axe du carottier)
        4 = bed_dip (dip - pendage de la stratification)
        5 = bed_dip_strike (str_ - direction de la stratification)

    Toutes les autres metadonnees (site/formation/age/lithologie...)
    sont ecrites "n.d" - a completer ensuite (STARpaleomag_Py, "Complete
    sample information...").

    BUG CORRIGE ici (present dans la version initiale du collegue) : un
    specimen ATRM apparait normalement PLUSIEURS FOIS dans un .ANI reel
    (une ligne par variante jackknife A0/A+/A-/A1/B1/A2/B2.../B6, 15 au
    total - verifie sur Chili_Briques_pmag.ANI : "11CA0101A" y apparait
    15 fois, toujours avec les MEMES cin/caz/dip/str) - sans
    deduplication, le .prmag resultant aurait contenu 15 blocs
    identiques pour ce meme specimen au lieu d'un seul. Seule la
    PREMIERE ligne rencontree pour chaque specimen est gardee (les
    valeurs d'orientation sont de toute facon identiques d'une variante
    a l'autre du meme specimen).

    Retourne le nombre de specimens (uniques) ecrits."""
    records = []
    seen_specimens = set()

    with open(ani_path, "r", encoding="iso-8859-1", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue

            fields = line.split()
            if len(fields) < 6 or fields[0].upper() != "D":
                continue

            specimen = fields[1]
            if specimen in seen_specimens:
                continue

            try:
                dip = float(fields[2])
                azimuth = float(fields[3])
                bed_dip = float(fields[4])
                bed_dip_strike = float(fields[5])
            except ValueError:
                continue

            seen_specimens.add(specimen)
            sample = _legacy_ani_specimen_to_sample(specimen)
            site = _legacy_ani_sample_to_site(sample)
            records.append((specimen, sample, site, azimuth, dip, bed_dip_strike, bed_dip))

    with open(prmag_path, "w", encoding="utf-8", newline="\n") as f:
        for specimen, sample, site, azimuth, dip, bed_dip_strike, bed_dip in records:
            f.write(
                f"specimen: {specimen}\tsample: {sample}\tsite: {site}\t"
                f"volume: {volume:.2f}\tmass: n.d\tlat: n.d\tlon: n.d\t"
                "elevation: n.d\tstratigraphic_height: n.d\tcomment: n.d\n"
            )
            f.write(
                f"azimuth: {azimuth}\tdip: {dip}\tdate: n.d\t"
                "magnetic_azimuth: n.d\tsolar_azimuth: n.d\torient_tool: n.d\n"
            )
            f.write(f"bed_dip_strike: {bed_dip_strike}\tbed_dip: {bed_dip}\n")
            f.write(
                "formation: n.d\tage: n.d\tgeologic_classes: n.d\t"
                "geologic_types: n.d\tlithologies: n.d\tlocation: n.d\t"
                "obs: n.d\tmethod_codes: n.d\n"
            )
            f.write(
                "step\tcod1\tcod2\tx\ty\tz\terror\tquality\tinstrument\ts\t"
                "treat_temp\ttreat_ac_field\ttreat_dc_strongfield\ttreat_dc_lowfield\t"
                "treat_dc_field_phi\ttreat_dc_field_theta\t"
                "method_codes\tinstrument_codes\ttreat_step_num\n"
            )
            f.write("\n")  # separateur de bloc (voir read_prmag_specimens)

    return len(records)
