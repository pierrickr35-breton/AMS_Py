"""
Lecture d'une contribution MagIC dans AMS_Py - demande explicite utilisateur
("est-ce possible de permettre la lecture des contributions MAGIC dans
AMS_Py"). Une contribution (fichier combine telecharge de MagIC, ou table
specimens.txt/measurements.txt avec ses voisines) -> .prmag (mesures VIDES,
meme principe que create_prmag_from_legacy_ani/_from_asc) + .pmagani.

Tenseurs :
- specimens.aniso_s ET measurements.aniso_s (le modele MagIC 3.0 prevoit les
  deux) ; un meme tenseur present dans les deux n'est garde qu'une fois.
- aniso_s est normalise a trace 1 dans MagIC, a trace 3 dans .pmagani.
- repere : aniso_tilt_correction -1 garde tel quel. Sans ligne -1 pour un
  specimen, le tenseur 0 (ou 100) est ramene en repere specimen par les
  rotations inverses de corfor_tensor/corpen_tensor : .pmagani est toujours
  en repere specimen (voir ams_selection.apply_orientation). Sans orientation
  disponible pour cela, le tenseur est ecarte et signale.
- code2 : N0/A0/AA/AI selon le type (aniso_type puis method_codes). Pour l'AMS,
  l'etat du tenseur vient des codes LT- : XF/YF/ZF (AF sur un axe), FX/FY/FZ
  (AF 3 axes, dernier axe), TF (tumbler), IZ/Z- (IRM selon +Z/-Z), RM (autre
  IRM), TH (thermique), F1 (AMS selon le champ alternatif), DI (tenseur
  difference, description "difference"), etape en mT, degC ou A/m.
Orientation (inverse de magic_export cote STARpaleomag_Py) : caz = azimuth+90,
cin = -dip, str_ = bed_dip_direction-90. Echantillon sans orientation :
caz=90, cin=0 (rotation neutre, voir corfor_tensor), signale.
sites.aniso_v1/v2/v3 -> section "#site mean tensor results".
"""

import csv
import math
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from ams_selection import (
    AMSMeasurement,
    cart,
    _format_pmagani_specimen_line,
    _insert_pmagani_line,
    _rota2,
    create_empty_pmagani_if_missing,
)
from ams_stats import shape_params

csv.field_size_limit(10 ** 8)

MU0 = 4e-7 * math.pi
_SEPARATOR = ">>>>>>>>>>"
_SIDE_TABLES = ("samples", "sites", "specimens", "measurements", "locations")
_GEOLOGY = {"formation": "formation", "age": "age", "age_low": "age_low", "age_high": "age_high",
            "age_unit": "age_unit", "geologic_classes": "geologic_classes",
            "geologic_types": "geologic_types", "lithologies": "lithologies", "location": "location"}


@dataclass
class MagicAMSImport:
    measurements: List[AMSMeasurement] = field(default_factory=list)
    prmag_records: List[dict] = field(default_factory=list)
    mean_lines: List[str] = field(default_factory=list)
    report: List[str] = field(default_factory=list)


def read_magic_tables(path: str) -> Dict[str, List[Dict[str, str]]]:
    """Fichier combine (tables separees par >>>>>>>>>>, BOM possible) ; une
    table seule (specimens.txt) est completee par ses voisines du dossier."""
    tables = _parse_tables(path)
    if len(tables) == 1:
        folder = os.path.dirname(path)
        for name in _SIDE_TABLES:
            side = os.path.join(folder, f"{name}.txt")
            if name not in tables and os.path.exists(side):
                tables.update(_parse_tables(side))
    return tables


def _parse_tables(path: str) -> Dict[str, List[Dict[str, str]]]:
    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
        text = f.read().replace("\r\n", "\n").replace("\r", "\n")
    tables: Dict[str, List[Dict[str, str]]] = {}
    for block in re.split(rf"^\s*{_SEPARATOR}.*$", text, flags=re.M):
        lines = [l for l in block.split("\n") if l.strip()]
        if len(lines) < 2 or not lines[0].lower().startswith("tab"):
            continue
        name = lines[0].split("\t")[-1].strip().lower()
        rows = csv.DictReader(lines[1:], delimiter="\t", quoting=csv.QUOTE_NONE)
        tables.setdefault(name, []).extend(
            {k.strip(): (v or "").strip() for k, v in r.items() if k} for r in rows)
    return tables


def _float(row: Dict[str, str], key: str) -> Optional[float]:
    try:
        v = float(row.get(key, ""))
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def _codes(row: Dict[str, str]) -> List[str]:
    return [c.strip() for c in row.get("method_codes", "").split(":") if c.strip()]


def _first_values(rows: List[Dict[str, str]], key: str) -> Dict[str, Dict[str, str]]:
    """Plusieurs lignes par nom (une par repere...) : premiere valeur non vide
    de chaque colonne."""
    out: Dict[str, Dict[str, str]] = {}
    for r in rows:
        name = r.get(key, "")
        if not name:
            continue
        merged = out.setdefault(name, {})
        for k, v in r.items():
            if v and not merged.get(k):
                merged[k] = v
    return out


def _main_code(row: Dict[str, str]) -> Optional[str]:
    """code2 du type d'anisotropie - None si ce n'est pas reconnu."""
    atype = row.get("aniso_type", "").upper()
    codes = _codes(row)
    if atype == "ATRM" or "LP-AN-TRM" in codes:
        return "A0"
    if atype == "AARM" or "LP-AN-ARM" in codes:
        return "AA"
    if atype == "AIRM" or "LP-AN-IRM" in codes:
        return "AI"
    if atype == "AMS" or "LP-AN-MS" in codes or "LP-X" in codes or "LP-X-H" in codes:
        return "N0"
    return None


def _state(row: Dict[str, str]) -> Tuple[str, int, Optional[float]]:
    """(code2, etape, champ A/m) d'un tenseur AMS, d'apres les codes LT-."""
    codes = _codes(row)
    if "difference" in row.get("description", "").lower():
        return "DI", 0, None
    ac = _float(row, "treat_ac_field")
    mt = int(round(ac * 1000.0)) if ac is not None else 0
    fac = _float(row, "meas_field_ac")
    if "LP-X-H" in codes and fac is not None:
        am = fac / MU0
        return "F1", int(round(am)), am
    if "LT-IRM" in codes:
        dc = _float(row, "treat_dc_field")
        theta = _float(row, "treat_dc_field_theta")
        step = int(round(dc * 1000.0)) if dc is not None else 0
        if theta == 90.0:
            return "IZ", step, None
        return ("Z-", step, None) if theta == -90.0 else ("RM", step, None)
    if "LT-AF-Z-TUMB" in codes:
        return "TF", mt, None
    for c in codes:
        m = re.fullmatch(r"LT-AF-Z-([XYZ])", c)
        if m:
            return f"{m.group(1)}F", mt, None
        m = re.fullmatch(r"LT-AF-Z-[XYZ]{2}([XYZ])", c)
        if m:
            return f"F{m.group(1)}", mt, None
    if "LT-AF-Z" in codes:
        return "FF", mt, None
    temp = _float(row, "treat_temp")
    if "LT-T-Z" in codes and temp is not None and temp > 300.0:
        return "TH", int(round(temp - 273.0)), None
    return "N0", 0, None


def _tensor_matrix(s: List[float]) -> np.ndarray:
    return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])


def _back_to_specimen(a: np.ndarray, tilt: str, cin: float, caz: float,
                      dip: float, strike: float) -> np.ndarray:
    """Inverse de apply_orientation : 100 -> 0 (corpen inverse) puis 0 -> -1
    (corfor inverse) ; les deux matrices sont orthogonales."""
    if tilt == "100":
        rot = _rota2(strike, 0.0, dip)
        a = rot.T @ a @ rot
    tm = np.zeros((3, 3))
    tm[:, 0] = cart(1.0, caz - 90.0, -cin)
    tm[:, 1] = cart(1.0, caz, 0.0)
    tm[:, 2] = cart(1.0, caz - 90.0, 90.0 - cin)
    return tm.T @ a @ tm


def _orientation(sample: Dict[str, str], site: Dict[str, str]) -> dict:
    az, dp = _float(sample, "azimuth"), _float(sample, "dip")
    bdd = _float(sample, "bed_dip_direction")
    bd = _float(sample, "bed_dip")
    if bd is None:
        bdd, bd = _float(site, "bed_dip_direction"), _float(site, "bed_dip")
    oriented = az is not None and dp is not None
    return dict(oriented=oriented, caz=(az + 90.0) % 360.0 if oriented else 90.0,
                cin=-dp if oriented else 0.0, dip=bd or 0.0,
                str_=((bdd - 90.0) % 360.0) if bdd is not None else 0.0)


def _site_mean_lines(sites: List[Dict[str, str]], n_by_site: Counter) -> Tuple[List[str], int]:
    lines, bad = [], 0
    for r in sites:
        if not r.get("aniso_v1"):
            continue
        code2 = _main_code(r) or "N0"
        axes = []
        for key in ("aniso_v1", "aniso_v2", "aniso_v3"):
            p = [x.strip() for x in r.get(key, "").split(":")]
            try:
                alpha = (float(p[6]), float(p[9])) if len(p) >= 10 else (0.0, 0.0)
                axes.append((float(p[0]), float(p[1]), float(p[2]), *alpha))
            except (ValueError, IndexError):
                axes = []
                break
        if len(axes) != 3:
            bad += 1
            continue
        taus = [a[0] for a in axes]
        scale = 3.0 / sum(taus) if sum(taus) else 1.0  # trace 1 (PmagPy) ou trace 3
        k1, k2, k3 = (t * scale for t in taus)
        sp = shape_params(k1, k2, k3)
        fields = [r.get("site", ""), code2, str(n_by_site[(r.get("site", ""), code2)]),
                  f"{k1:.6E}", f"{k2:.6E}", f"{k3:.6E}"]
        fields += [f"{a[i]:.2f}" for a in axes for i in (1, 2)]
        fields += [f"{a[i]:.3f}" for a in axes for i in (3, 4)]
        fields += [f"{sp[k]:.6g}" for k in ("P", "T", "L", "F", "Pprim")]
        tilt = r.get("aniso_tilt_correction", "")
        fields += ["1" if tilt == "-1" else tilt, "Y", '"imported from a MagIC contribution"']
        lines.append("\t".join(fields) + "\n")
    return lines, bad


def build_ams_import(path: str) -> MagicAMSImport:
    tables = read_magic_tables(path)
    out = MagicAMSImport()
    specimens = _first_values(tables.get("specimens", []), "specimen")
    samples = _first_values(tables.get("samples", []), "sample")
    sites = _first_values(tables.get("sites", []), "site")

    # Toutes les lignes avec un tenseur, measurements d'abord (observations)
    rows = []
    for table in ("measurements", "specimens"):
        rows += [(table, r) for r in tables.get(table, []) if r.get("aniso_s") and r.get("specimen")]
    if not rows:
        out.report.append("no aniso_s in the specimens or measurements tables")

    has_minus1 = {r["specimen"] for t, r in rows if r.get("aniso_tilt_correction", "-1") in ("-1", "")}
    seen, skipped = set(), Counter()
    kept_specimens: Dict[str, dict] = {}
    n_by_site: Counter = Counter()
    frames_back = 0
    unoriented = set()
    for table, r in rows:
        spec = r["specimen"]
        tilt = r.get("aniso_tilt_correction", "") or "-1"
        if spec in has_minus1 and tilt != "-1":
            continue  # meme tenseur en repere geographique/tectonique
        main = _main_code(r)
        if main is None:
            skipped["unknown anisotropy type"] += 1
            continue
        try:
            s = [float(x) for x in r["aniso_s"].split(":")]
            if len(s) != 6:
                raise ValueError
        except ValueError:
            skipped["unreadable aniso_s"] += 1
            continue
        spec_row = specimens.get(spec, {})
        sample_name = spec_row.get("sample") or r.get("sample") or re.sub(r"[A-Za-z]+$", "", spec)
        sample = samples.get(sample_name, {})
        site_name = sample.get("site") or spec_row.get("site") or sample_name
        site = sites.get(site_name, {})
        ori = _orientation(sample, site)
        a = _tensor_matrix(s)
        if tilt in ("0", "100"):
            if not ori["oriented"]:
                skipped[f"tilt {tilt} tensor without sample orientation"] += 1
                continue
            a = _back_to_specimen(a, tilt, ori["cin"], ori["caz"], ori["dip"], ori["str_"])
            frames_back += 1
        elif tilt != "-1":
            skipped[f"aniso_tilt_correction {tilt}"] += 1
            continue
        trace = float(np.trace(a))
        if not trace:
            skipped["zero trace"] += 1
            continue
        a = a * 3.0 / trace
        key = (spec, tuple(np.round([a[0, 0], a[1, 1], a[2, 2], a[0, 1], a[1, 2], a[0, 2]], 4)))
        if key in seen:
            continue  # meme tenseur dans measurements et specimens
        seen.add(key)
        # Dans specimens, les codes LT- decrivent le traitement de remanence du
        # specimen, pas l'etat de la mesure d'AMS : seul measurements les utilise.
        if main != "N0":
            code2, etape, fieldam = main, 0, None
        elif table == "measurements" or "difference" in r.get("description", "").lower():
            code2, etape, fieldam = _state(r)
        else:
            code2, etape, fieldam = "N0", 0, None
        k = _float(r, "susc_chi_volume")
        if k is None:
            k = _float(r, "aniso_s_mean")
        m = AMSMeasurement(
            id=spec, etape=etape, code2=code2, cin=ori["cin"], caz=ori["caz"], dip=ori["dip"], str_=ori["str_"],
            s=(k or 0.0) * 1.0e5,
            k11=a[0, 0], k22=a[1, 1], k33=a[2, 2], k12=a[0, 1], k23=a[1, 2], k13=a[0, 2],
            info=r.get("description", "").replace('"', "'").replace("\t", " ")[:200])
        m.field = fieldam
        m.sigma = _float(r, "aniso_s_sigma")
        m.ftest, m.ftest12, m.ftest23 = (_float(r, "aniso_ftest"), _float(r, "aniso_ftest12"),
                                         _float(r, "aniso_ftest23"))
        m.quality = r.get("aniso_ftest_quality") if r.get("aniso_ftest_quality") in ("g", "b") else None
        n = _float(r, "aniso_s_n_measurements")
        m.n_positions = int(n) if n is not None else None
        out.measurements.append(m)
        n_by_site[(site_name, main)] += 1
        if not ori["oriented"]:
            unoriented.add(sample_name)
        if spec not in kept_specimens:
            vol = _float(spec_row, "volume")
            geo = {k: site.get(col, "") for k, col in _GEOLOGY.items()}
            geo["lat"], geo["lon"] = site.get("lat", ""), site.get("lon", "")
            kept_specimens[spec] = dict(specimen=spec, sample=sample_name, site=site_name,
                                        volume=vol * 1.0e6 if vol else 10.8, **ori, **geo)

    out.prmag_records = list(kept_specimens.values())
    out.mean_lines, bad_means = _site_mean_lines(tables.get("sites", []), n_by_site)
    n_no_k = sum(1 for m in out.measurements if not m.s)
    if n_no_k:
        out.report.append(f"{n_no_k} tensor(s) without mean susceptibility in the contribution (s = 0)")
    if frames_back:
        out.report.append(f"{frames_back} tensor(s) brought back to specimen coordinates (no -1 row)")
    if unoriented:
        out.report.append(f"{len(unoriented)} sample(s) without azimuth/dip: neutral orientation, "
                          "in situ and tilt-corrected views are not meaningful for them")
    for why, n in skipped.items():
        out.report.append(f"{n} tensor(s) skipped: {why}")
    if bad_means:
        out.report.append(f"{bad_means} site mean(s) skipped: unreadable aniso_v1/v2/v3")
    n_sample_means = sum(1 for r in tables.get("samples", []) if r.get("aniso_v1"))
    if n_sample_means:
        out.report.append(f"{n_sample_means} sample mean(s) in the samples table not imported")
    return out


def _nd(v) -> str:
    return "n.d" if v in ("", None) else str(v)


def write_prmag(path: str, records: List[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for r in records:
            f.write(f"specimen: {r['specimen']}\tsample: {r['sample']}\tsite: {r['site']}\t"
                    f"volume: {r['volume']:.2f}\tmass: n.d\tlat: {_nd(r['lat'])}\tlon: {_nd(r['lon'])}\t"
                    "elevation: n.d\tstratigraphic_height: n.d\tcomment: n.d\n")
            f.write(f"azimuth: {r['caz']}\tdip: {r['cin']}\tdate: n.d\t"
                    "magnetic_azimuth: n.d\tsolar_azimuth: n.d\torient_tool: n.d\n")
            f.write(f"bed_dip_strike: {r['str_']}\tbed_dip: {r['dip']}\n")
            f.write(f"formation: {_nd(r['formation'])}\tage: {_nd(r['age'])}\tage_low: {_nd(r['age_low'])}\t"
                    f"age_high: {_nd(r['age_high'])}\tage_unit: {_nd(r['age_unit'])}\t"
                    f"geologic_classes: {_nd(r['geologic_classes'])}\tgeologic_types: {_nd(r['geologic_types'])}\t"
                    f"lithologies: {_nd(r['lithologies'])}\tlocation: {_nd(r['location'])}\t"
                    "obs: n.d\tmethod_codes: n.d\n")
            f.write("step\tcod1\tcod2\tx\ty\tz\terror\tquality\tinstrument\ts\t"
                    "treat_temp\ttreat_ac_field\ttreat_dc_strongfield\ttreat_dc_lowfield\t"
                    "treat_dc_field_phi\ttreat_dc_field_theta\t"
                    "method_codes\tinstrument_codes\ttreat_step_num\n\n")


def create_prmag_pmagani_from_magic(magic_path: str, prmag_path: str, pmagani_path: str) -> MagicAMSImport:
    """Ecrit les deux fichiers. Refuse d'ecraser un .pmagani existant (il peut
    contenir des moyennes ou un marquage d'export)."""
    if os.path.exists(pmagani_path):
        raise FileExistsError(pmagani_path)
    result = build_ams_import(magic_path)
    if not result.measurements:
        return result
    write_prmag(prmag_path, result.prmag_records)
    create_empty_pmagani_if_missing(pmagani_path)
    with open(pmagani_path, "a", encoding="utf-8", newline="\n") as f:
        for m in result.measurements:
            f.write(_format_pmagani_specimen_line(m))
    for line in result.mean_lines:
        _insert_pmagani_line(pmagani_path, line, is_mean=True)
    return result
