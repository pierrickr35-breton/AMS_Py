"""Import des fichiers .asc AGICO (rapport texte brut du kappabridge, logiciel
SUSAR/Anisoft) vers le format .pmagani - port de `readasc` (reference/
AMS_OSX_AWE/lect_asc.f:262-434) - demande explicite utilisateur ("can you also
change the import asc files to new format").

VERIFIE contre un vrai fichier .asc + son .ANI deja converti (Py_Dev/Caleu/
Caleu.ASC + Caleu.ASC.ANI, fournis par l'utilisateur) : 616/624 blocs
retrouves EXACTEMENT (id + tenseur + susceptibilite + cin/caz/dip/str, 0
divergence d'orientation) - les 8 restants sont soit un bloc reellement
corrompu dans le fichier source (deux specimens colles sans separateur de
page, ligne "Geograph  D ... 198PL0407B ..."), soit des specimens absents de
ce .ANI (mesures ajoutees plus tard - ex. "99PP1403B1"/"99PP1403B2" alors que
le .ANI ne connait qu'un "99PP1403B" avec un tenseur different).

DEUX ECARTS CONFIRMES entre la source Fortran fournie et le comportement REEL
observe (le .ANI reel ne correspond PAS a ce que lit_asc.f produirait tel
quel - cette source a manifestement diverge de la version qui a genere ce
fichier) :
  - `amsfichier(i).id` n'est JAMAIS assigne dans la source pour un specimen
    SANS suffixe _Kre/_KIm (seulement dans les 2 lignes conditionnelles) ;
    le vrai .ANI a pourtant le bon id partout -> ici, id=samplename TOUJOURS.
  - `amsfichier(i).etape=20` est inconditionnel dans la source ; le vrai .ANI
    a systematiquement etape=0 -> ici, etape=0 par defaut (le cas _KIm,
    ratio de susceptibilite, reste porte tel quel mais N'EST PAS verifie,
    aucun specimen _Kre/_KIm dans Caleu.ASC).

Branches NON exercees par Caleu.ASC (portees depuis le Fortran tel quel,
best-effort, PAS verifiees) : susceptibilite via marqueur F1/F3 (readings
multi-frequence), suffixes _Kre/_KIm, tenseur negatif (code2="I-").
"""

import os
from typing import List, Optional, Tuple

from ams_selection import AMSMeasurement, _PMAGANI_HEADER, _PMAGANI_UNITS_NOTE


def _parse_op_quad(line: str) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[int]]:
    seg = line[23:]
    try:
        return int(seg[0:2]), int(seg[4:6]), int(seg[8:10]), int(seg[12:14])
    except (ValueError, IndexError):
        return None, None, None, None


def parse_asc_file(path: str) -> Tuple[List[AMSMeasurement], List[str]]:
    """Retourne (measurements, warnings). Un bloc dont une section attendue
    (Azi/Dip/Specimen) est absente ou illisible est saute (warning ajoute),
    plutot que de faire echouer tout l'import - equivalent des `go to 14`/
    `go to 24` (abandon du bloc courant) du Fortran d'origine."""
    with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
        lines = [raw.rstrip("\n") for raw in f]

    out: List[AMSMeasurement] = []
    warnings: List[str] = []
    i, n = 0, len(lines)
    suscep0 = None  # derniere susceptibilite _Kre (pour le ratio _KIm)

    while i < n - 1:
        chaine, chaine2 = lines[i], lines[i + 1]
        if not chaine2.startswith("*"):
            i += 1
            continue

        block_start_line = i + 1  # 1-indexe, pour les messages
        raw = chaine[1:] if chaine.startswith("\x0c") else chaine
        samplename = raw[:14].strip()
        i += 2

        try:
            if not samplename or " " in samplename:
                # cas reel rencontre (Caleu.ASC) : un saut de page manquant
                # colle la fin d'un bloc au debut du suivant ("Geograph  D
                # ...198PL0407B...") - un id de specimen valide ne contient
                # jamais d'espace interne, donc un espace ici signale une
                # corruption du fichier source, pas une lecture correcte.
                raise ValueError(f"malformed specimen id {samplename!r} (corrupted block boundary?)")
            i += 1  # ligne vide apres '****'
            azi_line = lines[i]; i += 1
            if not azi_line.startswith("Azi"):
                raise ValueError("'Azi' line not found where expected")
            iaz = int(azi_line[6:9])
            ii1, ii2, ii3, ii4 = _parse_op_quad(azi_line)

            i += 1  # ligne vide
            dip_line = lines[i]; i += 1
            if not dip_line.startswith("Dip"):
                raise ValueError("'Dip' line not found where expected")
            icin = int(dip_line[6:9])

            cin = (90 - icin) if ii2 == 90 else float(icin)
            if ii3 == 6:
                caz = float(iaz - 90)
            elif ii3 == 12:
                caz = float(iaz + 90)
            else:  # ii3 == 3 (seul cas rencontre dans les donnees reelles)
                caz = float(iaz)

            strdip_line = None
            for _ in range(5):
                strdip_line = lines[i]; i += 1
            chars = list(strdip_line)
            for k in range(min(20, len(chars))):
                if chars[k] == "/":
                    chars[k] = " "
            toks = "".join(chars).split()
            if len(toks) < 3:
                raise ValueError("strike/dip line not parseable")
            istr, idip = int(toks[1]), int(toks[2])
            str_ = (istr - 90.0) if ii4 == 0 else float(istr)
            dip = float(idip)

            # Scan jusqu'au marqueur "  susc.  " (branche F1/F3 : non
            # verifiee, portee telle quelle depuis le Fortran)
            s_val, sigma_pct, ftest, ftest12, ftest23 = None, None, None, None, None
            truc2 = "N0"
            found_susc = False
            while i < n:
                line = lines[i]; i += 1
                if line[:9] == "  susc.  ":
                    i += 1  # ligne vide
                    data_line = lines[i]; i += 1
                    dtoks = data_line.split()
                    if not dtoks:
                        raise ValueError("susceptibility data line empty")
                    s_val = float(dtoks[0])
                    if len(dtoks) >= 5:
                        try:
                            sigma_pct = float(dtoks[2])
                            ftest = float(dtoks[3])
                            ftest12 = float(dtoks[4])
                            if len(dtoks) >= 6:
                                ftest23 = float(dtoks[5])
                        except ValueError:
                            pass
                    found_susc = True
                    break
                if len(line) > 6 and line[5:7] in ("F1", "F3"):
                    ftoks = line.split()
                    if len(ftoks) >= 3:
                        truc2 = ftoks[1]
                        s_val = float(ftoks[2])
                    found_susc = True
                    break
            if not found_susc or s_val is None:
                raise ValueError("susceptibility marker not found")
            s_val *= 100000.0

            for _ in range(6):
                i += 1
            conf_line = lines[i]; i += 1
            ctoks = conf_line.split()
            if len(ctoks) < 6:
                raise ValueError("confidence angles line not parseable")
            a1, a2, a3 = float(ctoks[3]), float(ctoks[4]), float(ctoks[5])

            k11 = k22 = k33 = k12 = k23 = k13 = None
            while i < n:
                line = lines[i]; i += 1
                if line[:8] == "Specimen":
                    k11 = float(line[45:52]); k22 = float(line[54:61]); k33 = float(line[63:70])
                    line2 = lines[i]; i += 1
                    k12 = float(line2[45:52]); k23 = float(line2[54:61]); k13 = float(line2[63:70])
                    break
            if k11 is None:
                raise ValueError("'Specimen' tensor line not found")

        except (ValueError, IndexError) as e:
            warnings.append(f"line {block_start_line}: skipped block ({e})")
            continue

        specimen_id = samplename
        etape = 0  # voir docstring du module : 0, pas 20 (verifie contre Caleu.ASC.ANI)
        nlen = len(samplename)
        if nlen > 4 and samplename[nlen - 4:] == "_Kre":
            truc2, specimen_id, suscep0 = "RE", samplename[:nlen - 4], s_val
        elif nlen > 4 and samplename[nlen - 4:] == "_KIm":
            truc2, specimen_id = "IM", samplename[:nlen - 4]
            if suscep0 and s_val:
                etape = min(9999, abs(int(suscep0 / s_val)))

        info = f"error angles: {a1:6.1f} {a2:6.1f} {a3:6.1f}"
        code2 = truc2
        if s_val < 0:
            info = f"negative K: {info.strip()}"
            code2 = "I-"
            k11, k22, k33, k12, k23, k13 = -k11, -k22, -k33, -k12, -k23, -k13

        out.append(AMSMeasurement(
            id=specimen_id, etape=etape, code2=code2.ljust(2)[:2],
            cin=cin, caz=caz, dip=dip, str_=str_,
            k11=k11, k22=k22, k33=k33, k12=k12, k23=k23, k13=k13, s=s_val,
            info=info, sigma=sigma_pct, ftest=ftest, ftest12=ftest12, ftest23=ftest23,
        ))

    return out, warnings


def import_asc_file(asc_path: str, new_path: Optional[str] = None) -> Tuple[str, List[AMSMeasurement], List[str]]:
    """Convertit un fichier .asc directement vers .pmagani (plus vers
    l'ancien .ANI) - demande explicite utilisateur ("change the import asc
    files to new format"). `new_path` par defaut : voir ams_prmag.
    ani_path_for(asc_path). cin/caz/dip/str_ SONT ecrits (ils viennent
    reellement du .asc, pas d'un .prmag) ; sigma/ftest/ftest12/ftest23 sont
    popules quand le fichier .asc les fournit (section "Tests for
    anisotropy") plutot que laisses "n.d"."""
    from ams_prmag import ani_path_for
    if new_path is None:
        new_path = ani_path_for(asc_path)
    measurements, warnings = parse_asc_file(asc_path)

    def fmt_stat(v):
        return "n.d" if v is None else f"{v:.6g}"

    with open(new_path, "w", encoding="utf-8") as out:
        out.write(f"# pmagani v1 - imported from ASC: {os.path.basename(asc_path)}\n")
        out.write(_PMAGANI_UNITS_NOTE)
        out.write("\t".join(_PMAGANI_HEADER) + "\n")
        for m in measurements:
            info = f'"{m.info}"' if m.info else '""'
            fields = [
                m.id, m.code2, str(m.etape),
                f"{m.k11:.6E}", f"{m.k22:.6E}", f"{m.k33:.6E}",
                f"{m.k12:.6E}", f"{m.k23:.6E}", f"{m.k13:.6E}",
                f"{m.s:.5f}", "n.d",
                fmt_stat(m.sigma), fmt_stat(m.ftest), fmt_stat(m.ftest12), fmt_stat(m.ftest23),
                "n.d",
                info,
            ]
            out.write("\t".join(fields) + "\n")
    return new_path, measurements, warnings
