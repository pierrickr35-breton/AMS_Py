"""
Diagrammes "Anisotropy parameters" (Flinn, P'/T) : contrairement au
stereonet (`ams_stereo.py`, port pixel-pres via `plotlib.PlotContext`), ces
graphiques sont reconstruits directement avec matplotlib (memes echelles/
graduations automatiques, meme choix que STARpaleomag_Py pour xygraph.py/
susceptibility.py) plutot que de reporter la logique d'ajustement iteratif
des graduations du Fortran (`flinn`/`flinax`/`tpprim`/`tpprimax`,
anisotropie.f) - dont le rapport d'exploration a confirme un bug reel de
troncature entiere sur les etiquettes d'axe (`fpn=1.+(i/10.)` avec division
entiere Fortran) qui n'a aucune raison d'etre reproduit ici.

Formules (L=k1/k2, F=k2/k3, T, P') identiques a `ams_stats.shape_params`,
deja verifiees exactes contre magicams.txt reel.

Les deux diagrammes sont empiles dans UNE seule figure (meme motif que
`xygraph.build_xygraph_figure` cote STARpaleomag_Py, panneaux `add_subplot(211)`/
`(212)`), un seul menu "Anisotropy parameters" (remplace les anciennes
entrees separees Flinn/T-Pprim). Chaque point a sa propre couleur (cycle
matplotlib par defaut) et apparait dans la legende (`label=id`, meme motif
que xygraph.py) - uniquement si moins de 15 entrees (sinon illisible).

Source des points : mesures brutes (`self.selection`) ET/OU resultats de
moyenne tensorielle deja calcules (`self.mean_results`, equivalent du choix
"donnees(0) resultats(1) d+r(2)" de `paramster`/`flinn`/`tpprim`) - voir
`shape_entries_from_measurements`/`shape_entries_from_mean_results`."""

from typing import List, Optional, Tuple

from matplotlib.figure import Figure

from ams_selection import AMSMeasurement
from ams_stats import principal_axes, shape_params, TensorialMeanResult

_LABEL_THRESHOLD = 15

ShapeEntry = Tuple[str, dict]  # (identifiant, dict de ams_stats.shape_params)


def shape_entries_from_measurements(
    measurements: List[AMSMeasurement], orientation: int,
) -> List[ShapeEntry]:
    """Equivalent du cas "donnees(0)" : un point par mesure brute de
    `measurements`, decomposee dans l'orientation courante."""
    out = []
    for m in measurements:
        axes = principal_axes(m, orientation)
        k1, k2, k3 = axes[0][0], axes[1][0], axes[2][0]
        out.append((m.id, shape_params(k1, k2, k3)))
    return out


def shape_entries_from_mean_results(results: List[TensorialMeanResult]) -> List[ShapeEntry]:
    """Equivalent du cas "resultats(1)" : un point par tenseur moyen deja
    calcule (`Calcul > Tensorial mean`, conserve dans `self.mean_results`)."""
    out = []
    for i, res in enumerate(results, start=1):
        if not res.axes:
            continue
        k1, k2, k3 = res.axes[0].eigenvalue, res.axes[1].eigenvalue, res.axes[2].eigenvalue
        label = res.id or f"mean #{i} (n={res.n})"
        out.append((label, shape_params(k1, k2, k3)))
    return out


def _draw_flinn(ax, entries: List[ShapeEntry]) -> list:
    """Retourne les `handles` (un par entree, meme ordre) - PAS de legende
    ici : une seule legende partagee (Flinn + P'/T) est construite par
    `build_anisotropy_parameters_figure`, ancree a droite de CE panneau."""
    handles = []
    if entries:
        fs = [sp["F"] for _id, sp in entries]
        ls = [sp["L"] for _id, sp in entries]
        for (eid, _sp), f, l in zip(entries, fs, ls):
            handles.append(ax.scatter([f], [l], marker="o", s=40, label=eid))
        lo = min(min(fs), min(ls), 1.0) - 0.01
        hi = max(max(fs), max(ls), 1.0) + 0.01
        ax.plot([lo, hi], [lo, hi], color="gray", linestyle="--", linewidth=0.7)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("F = k2/k3  (foliation)")
    ax.set_ylabel("L = k1/k2  (lineation)")
    ax.set_title("Flinn diagram")
    return handles


def _draw_pprim_t(ax, entries: List[ShapeEntry]) -> None:
    """P' en abscisse, T en ordonnee (P'/T, PAS T/P') - pas de legende
    (partagee avec Flinn, voir `_draw_flinn`), memes couleurs par entree
    (cycle matplotlib, meme ordre que Flinn -> couleurs coherentes entre
    les deux panneaux)."""
    if entries:
        ps = [sp["Pprim"] for _id, sp in entries]
        ts = [sp["T"] for _id, sp in entries]
        for (eid, _sp), p, t in zip(entries, ps, ts):
            ax.scatter([p], [t], marker="o", s=40, label=eid)
        ax.set_xlim(left=1.0)
    ax.axhline(0.0, color="gray", linestyle="--", linewidth=0.7)
    ax.set_ylim(-1.0, 1.0)
    ax.set_ylabel("T  (prolate ← 0 → oblate)")
    ax.set_xlabel("P'  (degree of anisotropy)")
    ax.set_title("P' / T diagram")


def build_anisotropy_parameters_figure(
    entries: List[ShapeEntry],
    fig: Optional[Figure] = None,
) -> Figure:
    """Equivalent (matplotlib, pas un port pixel-pres) de `flinn`+`tpprim`
    (menu Graphics > Anisotropy parameters, entree UNIQUE qui remplace les
    deux anciennes) reunis sur UNE page : Flinn (L=k1/k2 en ordonnee,
    F=k2/k3 en abscisse - meme orientation d'axes que le Fortran, `call
    symbol(y,x,...)`) en haut, P'/T (P' en abscisse, T en ordonnee) en bas.
    UNE SEULE legende des numeros d'echantillon, partagee par les deux
    graphiques, ancree a droite du diagramme de Flinn - pas une legende par
    panneau (demande explicite de l'utilisateur). Affichee seulement si
    moins de 15 entrees (sinon illisible)."""
    show_labels = len(entries) < _LABEL_THRESHOLD

    if fig is None:
        fig = Figure(figsize=(7.0, 9.5), dpi=100)
    else:
        fig.clear()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    handles = _draw_flinn(ax1, entries)
    _draw_pprim_t(ax2, entries)

    fig.suptitle(f"{len(entries)} donnee(s)" if entries else "")
    fig.tight_layout(rect=[0, 0, 0.78, 0.96] if (entries and show_labels) else None)
    if entries and show_labels:
        labels = [eid for eid, _sp in entries]
        ax1.legend(
            handles, labels, fontsize=7, loc="upper left",
            bbox_to_anchor=(1.03, 1.0), borderaxespad=0.0,
        )
    return fig
