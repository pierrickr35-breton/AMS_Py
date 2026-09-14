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

import math
from typing import List, Optional, Tuple

from matplotlib.figure import Figure

from ams_selection import AMSMeasurement, is_imaginary_component
from ams_stats import principal_axes, shape_params, TensorialMeanResult

_LABEL_THRESHOLD = 15

# (identifiant, dict de ams_stats.shape_params, vient du canal imaginaire ?)
# le 3e champ pilote le CHOIX DE PANNEAU (Im/Re separes, voir
# build_anisotropy_parameters_figure) - PAS le style du marqueur (le fill
# gris distinctif a ete retire des graphiques XY, reserve au stereonet -
# demande explicite utilisateur "in the various XY plots, it is not
# necessary to fill the symbols in grey. It is usefull only on the
# stereo plots to compare Re and Im"). Toujours False pour un tenseur
# moyen (pas de notion Re/Im au niveau d'une moyenne - voir
# shape_entries_from_mean_results).
ShapeEntry = Tuple[str, dict, bool]


def shape_entries_from_measurements(
    measurements: List[AMSMeasurement], orientation: int,
) -> List[ShapeEntry]:
    """Equivalent du cas "donnees(0)" : un point par mesure brute de
    `measurements`, decomposee dans l'orientation courante."""
    out = []
    for m in measurements:
        axes = principal_axes(m, orientation)
        k1, k2, k3 = axes[0][0], axes[1][0], axes[2][0]
        out.append((m.id, shape_params(k1, k2, k3), is_imaginary_component(m)))
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
        out.append((label, shape_params(k1, k2, k3), False))
    return out


# P=k1/k3 (degre d'anisotropie) au-dela duquel un point est considere
# comme un artefact/outlier (tenseur quasi-isotrope avec un axe bruite,
# specimen mal mesure...) plutot qu'une vraie anisotropie forte - demande
# explicite utilisateur ("consider as outliers samples with P > 5") : un
# seul point a P tres eleve ecrase sinon l'echelle des deux panneaux et
# rend illisibles tous les autres points.
P_OUTLIER_THRESHOLD = 5.0

# Attributs poses sur un artiste matplotlib (scatter d'un point, ou
# XAxis d'un panneau) pour le clic interactif cote app.py (_on_plot_pick,
# gestionnaire unique de pick_event partage par tous les graphiques -
# meme mecanisme que STARpaleomag_Py xygraph._pick_points/plotlib.
# PlotContext.pick_point, ce module restant volontairement hors
# PlotContext). "aniso_specimen" (demande explicite utilisateur - "click
# on point to have the specimen number") : `_ams_pick_data` porte l'id.
# "aniso_xaxis" (demande explicite utilisateur - "click on the X axis to
# change the max value") : pose sur l'artiste XAxis lui-meme (`set_
# picker` fait declencher pick_event des qu'un clic tombe pres de la
# ligne/graduations/etiquette de l'axe, pas seulement dans les donnees).


def _register_point_pick(h, specimen_id: str) -> None:
    h.set_picker(True)
    h._ams_pick_kind = "aniso_specimen"
    h._ams_pick_data = specimen_id


def _register_xaxis_pick(ax) -> None:
    ax.xaxis.set_picker(True)
    ax.xaxis._ams_pick_kind = "aniso_xaxis"
    ax.xaxis._ams_pick_data = None


def _draw_flinn(ax, entries: List[ShapeEntry], max_scale: Optional[float] = None) -> list:
    """Retourne les `handles` (un par entree, meme ordre) - PAS de legende
    ici : une seule legende partagee (Flinn + P'/T) est construite par
    `build_anisotropy_parameters_figure`, ancree a droite de CE panneau.
    Marqueur plein normal (fill=stroke=couleur du cycle) pour TOUS les
    points, Re comme Im - le fill gris distinctif reste reserve au
    stereonet (voir ams_stereo.draw_stereo_axes) : demande explicite
    utilisateur ("in the various XY plots, it is not necessary to fill
    the symbols in grey. It is usefull only on the stereo plots to
    compare Re and Im").

    `max_scale` (defaut None = auto sur les donnees, comme avant) : borne
    haute explicite choisie par l'utilisateur ("est-ce possible de
    selectionner l'echelle") pour F et L, le diagramme restant carre
    (memes bornes en x et y, ancre sur la droite 1:1). Chaque point est
    cliquable (id specimen, voir _register_point_pick) et l'axe X aussi
    (voir _register_xaxis_pick)."""
    handles = []
    if entries:
        fs = [sp["F"] for _id, sp, _im in entries]
        ls = [sp["L"] for _id, sp, _im in entries]
        for (eid, _sp, _im), f, l in zip(entries, fs, ls):
            h = ax.scatter([f], [l], marker="o", s=40, label=eid)
            _register_point_pick(h, eid)
            handles.append(h)
        if max_scale is not None:
            lo, hi = 1.0 - 0.01, max_scale
        else:
            lo = min(min(fs), min(ls), 1.0) - 0.01
            hi = max(max(fs), max(ls), 1.0) + 0.01
        ax.plot([lo, hi], [lo, hi], color="gray", linestyle="--", linewidth=0.7)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("F = k2/k3  (foliation)")
    ax.set_ylabel("L = k1/k2  (lineation)")
    ax.set_title("Flinn diagram")
    _register_xaxis_pick(ax)
    return handles


def _draw_pprim_t(
    ax, entries: List[ShapeEntry], max_scale: Optional[float] = None, title_suffix: str = "",
) -> list:
    """P' en abscisse, T en ordonnee (P'/T, PAS T/P'). Retourne les
    `handles` (un par entree, meme ordre) - utilise pour la legende
    partagee quand ce panneau n'a pas de Flinn a cote (voir
    `build_anisotropy_parameters_figure`, cas Re/Im separes). `max_scale`
    : meme borne haute explicite que `_draw_flinn`, appliquee ici a P'
    (T reste toujours -1..1, borne par definition). `title_suffix`
    distingue les deux panneaux P'/T quand ils sont empiles Im/Re (voir
    build_anisotropy_parameters_figure) plutot que Flinn+P'/T. Chaque
    point est cliquable (id specimen) et l'axe X aussi (changer
    l'echelle) - memes mecanismes que `_draw_flinn`. Marqueur plein
    normal pour tous les points (voir docstring `_draw_flinn` - fill gris
    retire, reserve au stereonet)."""
    handles = []
    if entries:
        ps = [sp["Pprim"] for _id, sp, _im in entries]
        ts = [sp["T"] for _id, sp, _im in entries]
        for (eid, _sp, _im), p, t in zip(entries, ps, ts):
            h = ax.scatter([p], [t], marker="o", s=40, label=eid)
            _register_point_pick(h, eid)
            handles.append(h)
        if max_scale is not None:
            ax.set_xlim(1.0, max_scale)
        else:
            ax.set_xlim(left=1.0)
    ax.axhline(0.0, color="gray", linestyle="--", linewidth=0.7)
    ax.set_ylim(-1.0, 1.0)
    ax.set_ylabel("T  (prolate ← 0 → oblate)")
    ax.set_xlabel("P'  (degree of anisotropy)")
    ax.set_title("P' / T diagram" + title_suffix)
    _register_xaxis_pick(ax)
    return handles


def build_anisotropy_parameters_figure(
    entries: List[ShapeEntry],
    fig: Optional[Figure] = None,
    max_scale_top: Optional[float] = None,
    max_scale_bottom: Optional[float] = None,
) -> Figure:
    """Equivalent (matplotlib, pas un port pixel-pres) de `flinn`+`tpprim`
    (menu Graphics > Anisotropy parameters, entree UNIQUE qui remplace les
    deux anciennes) reunis sur UNE page : Flinn (L=k1/k2 en ordonnee,
    F=k2/k3 en abscisse - meme orientation d'axes que le Fortran, `call
    symbol(y,x,...)`) en haut, P'/T (P' en abscisse, T en ordonnee) en bas.

    Exception (demande explicite utilisateur - "when we have samples in
    imaginary and real, instead of plotting the Flinn graph, plot P'T of
    real below and P'T of imaginary above") : des que `entries` melange
    au moins un point reel et un point imaginaire (voir `ShapeEntry`), le
    panneau Flinn est remplace par un DEUXIEME P'/T - Imaginary en haut,
    Real en bas (memes echelles/mise a l'echelle des outliers que le cas
    normal). Un jeu homogene (tout Re, tout Im, ou des moyennes - jamais
    marquees imaginaires, voir shape_entries_from_mean_results) garde la
    disposition Flinn+P'/T habituelle.

    UNE SEULE legende des numeros d'echantillon, partagee par les deux
    graphiques, ancree a droite du panneau du haut - pas une legende par
    panneau (demande explicite de l'utilisateur). Affichee seulement si
    moins de 15 entrees (sinon illisible).

    Les points avec P=k1/k3 > `P_OUTLIER_THRESHOLD` sont exclus des DEUX
    panneaux (et de leur mise a l'echelle automatique) - demande explicite
    utilisateur ("consider as outliers samples with P > 5") ; leur nombre
    est rapporte dans le titre plutot que silencieusement disparu.

    `max_scale_top`/`max_scale_bottom` : bornes hautes explicites
    INDEPENDANTES pour le panneau du haut et celui du bas - demande
    explicite utilisateur ("select the scale", puis correction "the
    scale should not apply to both plots (Re et Im) as they have very
    different ranges" : Im et Re ont des degres d'anisotropie de nature
    differente, une echelle commune ecrase systematiquement l'un des
    deux). Chacune peut aussi etre changee en cliquant sur l'axe X du
    panneau correspondant (voir _register_xaxis_pick / app._on_plot_pick,
    qui identifie le panneau par sa position haut/bas dans la figure,
    pas par son contenu - valable aussi bien pour Flinn/P'-T que pour
    P'-Im/P'-Re)."""
    kept = [e for e in entries if not (e[1]["P"] > P_OUTLIER_THRESHOLD)]
    n_outliers = len(entries) - len(kept)
    show_labels = len(kept) < _LABEL_THRESHOLD

    if fig is None:
        fig = Figure(figsize=(7.0, 9.5), dpi=100)
    else:
        fig.clear()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    im_entries = [e for e in kept if e[2]]
    re_entries = [e for e in kept if not e[2]]
    split_re_im = bool(im_entries) and bool(re_entries)

    if split_re_im:
        handles_top = _draw_pprim_t(ax1, im_entries, max_scale_top, title_suffix=" - Imaginary")
        handles_bottom = _draw_pprim_t(ax2, re_entries, max_scale_bottom, title_suffix=" - Real")
        handles = handles_top + handles_bottom
        legend_source = im_entries + re_entries
    else:
        handles = _draw_flinn(ax1, kept, max_scale_top)
        handles += _draw_pprim_t(ax2, kept, max_scale_bottom)
        legend_source = kept + kept

    title = f"{len(kept)} donnee(s)" if kept else ""
    if n_outliers:
        title += f"  ({n_outliers} outlier(s) P>{P_OUTLIER_THRESHOLD:g} excluded)"
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0, 0.78, 0.96] if (kept and show_labels) else None)
    if kept and show_labels:
        seen = set()
        leg_handles, leg_labels = [], []
        for (eid, _sp, _im), h in zip(legend_source, handles):
            if eid in seen:
                continue
            seen.add(eid)
            leg_handles.append(h)
            leg_labels.append(eid)
        ax1.legend(
            leg_handles, leg_labels, fontsize=7, loc="upper left",
            bbox_to_anchor=(1.03, 1.0), borderaxespad=0.0,
        )
    return fig


# ----------------------------------------------------------------------
# "Susceptibility / anisotropy degree" : susceptibilite (X, log) vs degre
# d'anisotropie P' (Y) - demande explicite utilisateur ("built a plot with
# two graphs of susceptibility on X and anisotropy degree on Y, above the
# Im and below the RE"). Un point par MESURE brute (pas de "resultats de
# moyenne" ici : une moyenne n'a pas de susceptibilite bulk individuelle
# comparable, voir TensorialMeanResult - ce plot n'a donc pas le choix
# "donnees(0) resultats(1) d+r(2)" de build_anisotropy_parameters_figure).
# Panneau du BAS = Real (susceptibilite toujours positive, log croissant
# normalement gauche->droite, "from min to max").
#
# Panneau du HAUT = Imaginary, EN DEUX MOITIES cote a cote partageant le
# meme axe Y (P') mais avec des echelles X INDEPENDANTES (demande
# explicite utilisateur, correction de la premiere version qui melangeait
# K+ et K- sur UN SEUL axe via |K| : "when K is positive, plot on the
# right panel and when negative to the left. We need two independant
# scales (and not with zero because of the log)" - une echelle log ne
# peut de toute facon pas traverser 0, donc pas de "vraie" jonction
# continue possible entre K+ et K- sur un meme axe) : moitie GAUCHE = K
# negatif, |K| en log, axe INVERSE (grandeur croissante VERS LA GAUCHE,
# en s'eloignant du centre - demande explicite anterieure conservee "for
# the negative K, plotting with the scale increasing to the left") ;
# moitie DROITE = K positif, K en log, axe normal (grandeur croissante
# vers la droite, en s'eloignant du centre) - meme filtrage outliers
# (P>5/indefini, voir P_OUTLIER_THRESHOLD) que build_anisotropy_
# parameters_figure, pour la meme raison (un point a degre d'anisotropie
# demesure ecrase l'echelle Y de tout le reste, PARTAGEE par les 2
# moities Imaginary).
# -----------------------------------------------------------------------

_SUSC_XLIM_MARGIN = 1.15  # facteur log au-dela du min/max exact - voir _draw_susceptibility_panel


def _susceptibility_degree_points(
    measurements: List[AMSMeasurement], orientation: int, degree_key: str = "Pprim",
) -> Tuple[List[Tuple[str, float, float]], int]:
    """Un point (id, susceptibilite SI SIGNEE, degre d'anisotropie) par
    mesure de `measurements`, filtrant les outliers (P>P_OUTLIER_
    THRESHOLD ou indefini, voir ams_stats.shape_params) - retourne aussi
    leur nombre. Le signe de la susceptibilite est CONSERVE (pas de
    valeur absolue ici) : c'est lui qui pilote le cote (gauche/droite) du
    panneau Imaginary dans build_susceptibility_anisotropy_figure - pour
    Real (toujours positif en pratique), le signe n'a simplement aucun
    effet. `degree_key` : cle de `shape_params` a utiliser comme
    "anisotropy degree" - "Pprim" (P', degre corrige de Jelinek) par
    defaut, meme convention que le panneau P'/T de build_anisotropy_
    parameters_figure ("P" brut reste un choix valide si demande
    explicitement)."""
    points = []
    n_outliers = 0
    for m in measurements:
        if m.s == 0.0:
            continue
        axes = principal_axes(m, orientation)
        sp = shape_params(axes[0][0], axes[1][0], axes[2][0])
        degree = sp[degree_key]
        if math.isnan(degree) or sp["P"] > P_OUTLIER_THRESHOLD:
            n_outliers += 1
            continue
        points.append((m.id, m.s * 1.0e-5, degree))
    return points, n_outliers


def _draw_susceptibility_panel(
    ax, points: List[Tuple[str, float, float]], title: str, reverse_x: bool,
    xlabel: str = "susceptibility (SI)", ylabel: str = "P'  (degree of anisotropy)",
) -> None:
    """Trace UN panneau/demi-panneau : `points` deja filtres au signe
    voulu et passes en |susceptibilite| (log ne supporte pas 0/negatif -
    voir build_susceptibility_anisotropy_figure pour le choix du cote)."""
    if points:
        xs = [x for _id, x, _y in points]
        for eid, x, y in points:
            h = ax.scatter([x], [y], marker="o", s=40, label=eid)
            _register_point_pick(h, eid)
        ax.set_xscale("log")
        # marge multiplicative (log) au-dela du min/max exact des points -
        # demande explicite utilisateur ("extend a little the min max
        # scales, sinon les symboles sont coupes en deux") : des points
        # au min/max EXACT se retrouvent avec leur marqueur coupe par le
        # cadre de l'axe. Meme facteur gere aussi, en passant, le cas
        # d'un seul point (lo==hi) - matplotlib "corrige" sinon des
        # bornes identiques en perdant le sens (`reverse_x`) au passage.
        lo, hi = min(xs) / _SUSC_XLIM_MARGIN, max(xs) * _SUSC_XLIM_MARGIN
        ax.set_xlim((hi, lo) if reverse_x else (lo, hi))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def build_susceptibility_anisotropy_figure(
    measurements: List[AMSMeasurement],
    orientation: int,
    fig: Optional[Figure] = None,
    degree_key: str = "Pprim",
) -> Figure:
    """"Susceptibility / anisotropy degree" : voir le commentaire de
    section ci-dessus. Separe `measurements` en Imaginary (haut) et Real
    (bas) via `is_imaginary_component` - un groupe vide laisse simplement
    son panneau nu (pas d'erreur), utile pour une selection homogene.
    Le panneau Imaginary est lui-meme scinde en 2 sous-panneaux cote a
    cote et JOINTIFS (K negatif a gauche, K positif a droite), memes
    bornes Y (`sharey`) mais echelles X INDEPENDANTES - une ligne
    verticale a leur jonction separe les deux groupes, graduations Y non
    dupliquees a droite - plutot qu'un seul axe continu, puisqu'un axe
    log ne peut de toute facon pas franchir 0 pour relier K- et K+."""
    im_meas = [m for m in measurements if is_imaginary_component(m)]
    re_meas = [m for m in measurements if not is_imaginary_component(m)]
    im_points, im_outliers = _susceptibility_degree_points(im_meas, orientation, degree_key)
    re_points, re_outliers = _susceptibility_degree_points(re_meas, orientation, degree_key)
    im_neg = [(eid, -x, y) for eid, x, y in im_points if x < 0.0]
    im_pos = [(eid, x, y) for eid, x, y in im_points if x > 0.0]

    if fig is None:
        fig = Figure(figsize=(7.0, 9.5), dpi=100)
    else:
        fig.clear()
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], wspace=0.0, hspace=0.35)
    ax_neg = fig.add_subplot(gs[0, 0])
    ax_pos = fig.add_subplot(gs[0, 1], sharey=ax_neg)
    ax_re = fig.add_subplot(gs[1, :])

    _draw_susceptibility_panel(
        ax_neg, im_neg, "Imaginary - negative K", reverse_x=True,
        xlabel="|susceptibility| (SI)")
    _draw_susceptibility_panel(
        ax_pos, im_pos, "Imaginary - positive K", reverse_x=False,
        xlabel="susceptibility (SI)", ylabel="")
    _draw_susceptibility_panel(ax_re, re_points, "Real", reverse_x=False)

    # jointure visuelle "axe casse" entre les 2 moities Imaginary : les 2
    # sous-panneaux se touchent (`wspace=0`), les graduations Y de droite
    # (dupliquees via sharey) sont masquees pour lire les 2 moities comme
    # UN panneau - et une ligne verticale distincte marque la frontiere
    # K-/K+ exactement a leur jonction (demande explicite utilisateur
    # "mettre une ligne verticale pour separer les groupes avec Im
    # negatives et positives") : seul le spine de GAUCHE est garde
    # visible (celui de droite serait rigoureusement au meme endroit -
    # le dedoubler ne ferait qu'epaissir le trait sans rien ajouter).
    ax_neg.spines["right"].set_visible(True)
    ax_neg.spines["right"].set_color("0.35")
    ax_neg.spines["right"].set_linewidth(1.3)
    ax_pos.spines["left"].set_visible(False)
    ax_pos.tick_params(left=False, labelleft=False)

    def _part(label, points, outliers):
        s = f"{label}: {len(points)}"
        if outliers:
            s += f" ({outliers} outlier(s) P>{P_OUTLIER_THRESHOLD:g} excluded)"
        return s

    fig.suptitle("  |  ".join(
        p for p in (
            _part("Im", im_points, im_outliers) if (im_points or im_outliers) else "",
            _part("Re", re_points, re_outliers) if (re_points or re_outliers) else "",
        ) if p
    ))
    # `fig.tight_layout()` emet un avertissement inoffensif mais bruyant
    # avec des Axes `sharey` issus d'un GridSpec (voir _draw_susceptibility_
    # panel/ax_pos ci-dessus) - marges fixes equivalentes via
    # `subplots_adjust`, sans cet avertissement.
    fig.subplots_adjust(left=0.11, right=0.97, top=0.93, bottom=0.07, wspace=0.05, hspace=0.35)
    return fig


# ----------------------------------------------------------------------
# "Im vs Re susceptibility" : Im susceptibilite (X) vs Re susceptibilite
# (Y), un point par SPECIMEN ayant a la fois une mesure Real et une
# mesure Imaginary - demande explicite utilisateur ("a final plot with
# the Im susceptibility on X and the Re susceptibility on Y same
# conventions as above for Kim + and -") : memes conventions que
# build_susceptibility_anisotropy_figure pour le cote Im (K negatif a
# gauche/K positif a droite, echelles log INDEPENDANTES, marge, ligne de
# separation) - mais Y (susceptibilite Re, PAS un ratio borne comme P')
# est ICI AUSSI en echelle log, partagee (`sharey`) entre les 2 moities.
# -----------------------------------------------------------------------

def _re_im_susceptibility_pairs(
    measurements: List[AMSMeasurement],
) -> List[Tuple[str, float, float]]:
    """Un point (id, susceptibilite Im SIGNEE, susceptibilite Re) par
    specimen present a la fois cote Real ET Imaginary de `measurements` -
    susceptibilite scalaire (`m.s`), INDEPENDANTE de l'orientation
    (contrairement aux axes propres/degre d'anisotropie), donc pas de
    parametre orientation ici. Un specimen n'ayant qu'un seul des deux
    canaux (mesure isolee, import partiel...) est silencieusement ignore
    - rien a comparer."""
    re_by_id, im_by_id = {}, {}
    for m in measurements:
        if m.s == 0.0:
            continue
        if is_imaginary_component(m):
            im_by_id.setdefault(m.id, m)
        else:
            re_by_id.setdefault(m.id, m)
    return [
        (id_, im_by_id[id_].s * 1.0e-5, re_m.s * 1.0e-5)
        for id_, re_m in re_by_id.items() if id_ in im_by_id
    ]


def _draw_re_im_panel(
    ax, points: List[Tuple[str, float, float]], title: str, reverse_x: bool, xlabel: str,
) -> list:
    handles = []
    if points:
        xs = [x for _id, x, _y in points]
        for eid, x, y in points:
            h = ax.scatter([x], [y], marker="o", s=40, label=eid)
            _register_point_pick(h, eid)
            handles.append(h)
        ax.set_xscale("log")
        lo, hi = min(xs) / _SUSC_XLIM_MARGIN, max(xs) * _SUSC_XLIM_MARGIN
        ax.set_xlim((hi, lo) if reverse_x else (lo, hi))
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    return handles


def build_im_re_susceptibility_figure(
    measurements: List[AMSMeasurement],
    fig: Optional[Figure] = None,
) -> Figure:
    """"Im vs Re susceptibility" : voir le commentaire de section
    ci-dessus. Scinde en 2 sous-panneaux JOINTIFS (K negatif a gauche,
    K positif a droite), memes bornes Y (`sharey`, susceptibilite Re en
    log) mais echelles X INDEPENDANTES, ligne verticale a leur jonction -
    mêmes conventions que build_susceptibility_anisotropy_figure."""
    pairs = _re_im_susceptibility_pairs(measurements)
    neg = [(eid, -x, y) for eid, x, y in pairs if x < 0.0]
    pos = [(eid, x, y) for eid, x, y in pairs if x > 0.0]
    show_labels = len(pairs) < _LABEL_THRESHOLD

    if fig is None:
        fig = Figure(figsize=(7.0, 5.5), dpi=100)
    else:
        fig.clear()
    gs = fig.add_gridspec(1, 2, wspace=0.0)
    ax_neg = fig.add_subplot(gs[0, 0])
    ax_pos = fig.add_subplot(gs[0, 1], sharey=ax_neg)

    handles_neg = _draw_re_im_panel(
        ax_neg, neg, "Im negative K", reverse_x=True, xlabel="|Im susceptibility| (SI)")
    handles_pos = _draw_re_im_panel(
        ax_pos, pos, "Im positive K", reverse_x=False, xlabel="Im susceptibility (SI)")
    ax_neg.set_ylabel("Re susceptibility (SI)")
    ax_pos.tick_params(left=False, labelleft=False)

    all_ys = [y for _id, _x, y in pairs]
    if all_ys:
        ax_neg.set_yscale("log")
        ax_neg.set_ylim(min(all_ys) / _SUSC_XLIM_MARGIN, max(all_ys) * _SUSC_XLIM_MARGIN)

    # meme ligne de separation K-/K+ que build_susceptibility_anisotropy_
    # figure (voir son commentaire pour le detail du choix).
    ax_neg.spines["right"].set_visible(True)
    ax_neg.spines["right"].set_color("0.35")
    ax_neg.spines["right"].set_linewidth(1.3)
    ax_pos.spines["left"].set_visible(False)

    fig.suptitle(f"{len(pairs)} specimen(s) with both Re and Im")
    if pairs and show_labels:
        handles = handles_neg + handles_pos
        labels = [eid for eid, _x, _y in (neg + pos)]
        ax_pos.legend(
            handles, labels, fontsize=7, loc="upper left",
            bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0,
        )
        fig.subplots_adjust(left=0.12, right=0.80, top=0.90, bottom=0.12, wspace=0.0)
    else:
        fig.subplots_adjust(left=0.12, right=0.97, top=0.90, bottom=0.12, wspace=0.0)
    return fig
