"""
Stereonet AMS : port de `stereo`/`stereo2`/`circle`/`stecor` (anisotropie.f)
via `plotlib.PlotContext` - meme bibliotheque graphique de bas niveau que
STARpaleomag_Py (confirme).

Contrairement a une premiere version qui empruntait la projection de
STARpaleomag_Py (`superc`, geree pour un pole de projection oblique), la lecture du
VRAI source de `circle` (anisotropie.f:2249-2335) montre que le trace du
reseau AMS passe par `stecor` (formule simple, SANS re-projection sur pole
oblique - `axcha`/`newcor`, non portes, restent hors perimetre) ET par
`symbo1` (divers.f:32-64) - PAS un appel direct a `plot`/`symbol`. C'est
`symbo1` qui pose probleme si on l'ignore : il (1) permute les coordonnees
(x,y)->(-y,x) et (2) REMAPPE le type de symbole demande vers un type
different reellement trace (table exacte ci-dessous, transcrite du
source). `stereo`/`stereo2` (trace des axes k1/k2/k3) passent par le MEME
`symbo1` - la table de remappage y est donc tout aussi determinante.

Simplification assumee (demande explicite de l'utilisateur, "evidemment
simplifier") : seul le pole de projection standard (equivalent la=-90 cote
STARpaleomag_Py) est supporte - le mode "projection non standard" de `paramster`
(pole personnalise via `axcha`/`newcor`) n'est pas porte."""

import math
from typing import List, Optional, Tuple

import numpy as np
from matplotlib.figure import Figure

from ams_selection import AMSMeasurement, is_imaginary_component
from ams_stats import principal_axes, TensorialMeanResult
from ams_bootstrap import BootstrapResult
from plotlib import PlotContext


def stecor(dec_deg: float, inc_deg: float, r: float, ams_iproj: int = 0) -> Tuple[float, float]:
    """Equivalent EXACT de `stecor` (anisotropie.f:2231-2241, formule
    verifiee sur le rapport d'exploration) : `ams_iproj` 1=stereographique,
    sinon (0, defaut) equiaire/Lambert - MEME convention que le prompt
    `paramster` ("projection stereo(1) lambert(0)"), pas celle de STARpaleomag_Py."""
    pi = math.pi
    d = math.radians(dec_deg)
    i = math.radians(inc_deg)
    if ams_iproj == 1:
        t = r * math.tan(pi / 4.0 - i / 2.0)
    else:
        t = r * math.sqrt(2.0) * math.sin(pi / 4.0 - i / 2.0)
    return t * math.cos(d), -t * math.sin(d)


# Table exacte de `symbo1` (divers.f:49-58) - remappe le type de symbole
# DEMANDE (cote de l'appelant, `stereo`/`circle`) vers celui REELLEMENT
# trace. Types plotlib (voir plotlib.py) : 8=cercle, 9=carre, 10=triangle,
# 14=cercle plein, 15=carre plein, 16=triangle plein, 5=plus.
_SYMBO1_MAP = {3: 5, 11: 15, 5: 9, 7: 14, 2: 16, 10: 8, 4: 14, 14: 8, 19: 8, 20: 10}


def _plott(ctx: PlotContext, x: float, y: float, mode: int, linestyle: str = "solid") -> None:
    """Equivalent de `plott` (divers.f:7-12) : `plot(-y,x,mode)`. `linestyle`
    (voir PlotContext.plot) transmis tel quel - defaut "solid", aucun
    effet sur les appelants qui ne le passent pas."""
    ctx.plot(-y, x, mode, linestyle=linestyle)


def _symbo1(ctx: PlotContext, x: float, y: float, h: float, ityp: int, nt: int = -1) -> None:
    """Equivalent simplifie de `symbo1` (divers.f:32-64) - PAS l'angle de
    rotation (`plotlib.PlotContext.symbol` ne le prend deja pas en charge,
    meme simplification que le reste du port)."""
    ib = _SYMBO1_MAP.get(ityp, ityp)
    ctx.symbol(-y, x, h, ib, nt)


_ORIENT_LABELS = {1: "Sample Coor.", 2: "In situ", 3: "Tilt corrected"}


def draw_stereo_net(
    ctx: PlotContext, orientation: int = 2, ams_iproj: int = 0,
    dimster: float = 15.0, show_orient_label: bool = True,
) -> float:
    """Equivalent de `circle` (anisotropie.f:2249-2335) : cercle du reseau
    (trace direct, invariant par rotation), graduations tous les 10 deg sur
    les 4 azimuts cardinaux + petits arcs 0-4 deg (`plott`), marqueur du
    pole (`symbo1`, type 3 -> "+"), N/E/S/W et etiquette d'orientation
    (`plottxt` DIRECT - PAS de wrapper, verifie identique a la version
    STARpaleomag_Py deja portee). Retourne le rayon `r` (cm)."""
    r = dimster / 3.0

    ctx.thickn(1.0)
    ctx.newpen(1)
    ctx.circl2(0.0, 0.0, r, 1, 0)

    ctx.thickn(0.25)
    for i in (90, 180, 270, 360):
        x, y = stecor(i, 0.0, r, ams_iproj)
        _plott(ctx, x, y, 3)
        x, y = stecor(i, 4.0, r, ams_iproj)
        _plott(ctx, x, y, 2)
        for j in range(0, 81, 10):
            x, y = stecor(float(i), float(j), r, ams_iproj)
            _symbo1(ctx, x, y, r / 30.0, 3)

    x, y = stecor(0.0, 90.0, r, ams_iproj)
    _symbo1(ctx, x, y, r / 7.0, 3)

    if show_orient_label:
        text = _ORIENT_LABELS.get(orientation, "")
        ctx.newpen(4)
        ctx.plottxt(-r - 4 * dimster / 84, (r + dimster / 21.6) - dimster / 27.63, dimster / 38.158692, text)
        ctx.newpen(1)

    ctx.plottxt(-dimster / 68.76, (r + dimster / 21.6) - dimster / 27.63, dimster / 38.158692, "N")
    ctx.plottxt(r + dimster / 190.9, -dimster / 65.6233, dimster / 38.158692, "E")
    ctx.plottxt(-r - 4 * dimster / 84, -dimster / 65.6233, dimster / 38.158692, "W")
    ctx.plottxt(-dimster / 68.76, (-r - dimster / 1050.0) - dimster / 27.63, dimster / 38.158692, "S")

    return r


# ipos/neg (anisotropie.f, `stereo`/`stereo2`) : types DEMANDES a symbo1 -
# une fois remappes (_SYMBO1_MAP), donnent carre(k1)/triangle(k2)/cercle(k3)
# pleins pour un axe "positif", ouverts pour un axe negatif non inverse.
_IPOS = {1: 11, 2: 2, 3: 7}   # -> plein : carre(15) / triangle(16) / cercle(14)
_NEG = {1: 5, 2: 20, 3: 19}   # -> ouvert : carre(9) / triangle(10) / cercle(8)
# BUG REEL corrige ici (demande explicite utilisateur, "is it possible to
# keep the Kmin in green as in the old app") : verifie directement contre
# le VRAI source Fortran (reference/AMS_OSX_AWE/anisotropie.f, subroutine
# `stereo` ligne ~2097-2103 ET `tratm` ligne ~2431-2433) - dans les DEUX,
# `ipen=i*2+1` PUIS un override explicite `if(i==3) ipen=4` : i=1(kmax)
# -> pen 3 (rouge), i=2(kint) -> pen 5 (bleu, PAS d'override), i=3(kmin)
# -> pen 4 (VERT, via l'override). Le port precedent utilisait la formule
# SANS l'override (pen 4 pour kint, pen 5 pour kmin) - kint et kmin
# avaient donc leurs couleurs INTERVERTIES par rapport a l'appli
# d'origine depuis le debut du port.
_AXIS_PEN = {1: 3, 2: 5, 3: 4}  # rouge(kmax)/bleu(kint)/vert(kmin)
_AXIS_COLOR = {1: "red", 2: "blue", 3: "green"}
# RGB exacts des memes couleurs (matplotlib "red"/"blue"/"green"), pour
# newpencol : point Im = stroke couleur d'axe + fill gris (au lieu de
# stroke=fill=couleur d'axe pour un point Re) - demande explicite
# utilisateur ("mettre les tenseurs d'Im avec un fill en gris et le
# stroke color rouge vert bleu ... pour differencier les Re des Im").
_AXIS_RGB = {1: (255, 0, 0), 2: (0, 0, 255), 3: (0, 128, 0)}
_IM_FILL_RGB = (160, 160, 160)
# ratio symbole moyenne/mesure individuelle - voir draw_stereo_mean_results
_MEAN_SIZE_RATIO = 5.0 / 3.0


def draw_stereo_axes(
    ctx: PlotContext,
    measurements: List[AMSMeasurement],
    orientation: int,
    r: float,
    invert_negative: bool = True,
    ams_iproj: int = 0,
    point_size: float = 0.324,
) -> None:
    """Equivalent du trace des axes k1/k2/k3 dans `stereo`/`stereo2`
    (anisotropie.f:~2090-2165) : pour chaque mesure, projette les 3 axes
    propres (NON pre-retournes - `principal_axes(..., flip_negative_
    inclination=False)`) via `stecor`+`symbo1`. `invert_negative` = 'iinv'
    du Fortran : si actif, un axe d'inclinaison negative est retourne a
    l'antipode (meme symbole "positif") ; sinon il reste a sa position
    reelle avec le symbole "neg" (ouvert)."""
    ctx.thickn(0.5)
    for m in measurements:
        im = is_imaginary_component(m)
        axes = principal_axes(m, orientation, flip_negative_inclination=False)
        for j, (_ev, dec, inc) in enumerate(axes, start=1):
            sym = _IPOS[j]
            if inc < 0.0:
                inc = -inc
                if invert_negative:
                    dec += 180.0
                    if dec > 360.0:
                        dec -= 360.0
                else:
                    sym = _NEG[j]
            x, y = stecor(dec, inc, r, ams_iproj)
            if im:
                sr, sg, sb = _AXIS_RGB[j]
                ctx.newpencol(sr, sg, sb, *_IM_FILL_RGB)
            else:
                ctx.newpen(_AXIS_PEN[j])
            _symbo1(ctx, x, y, point_size, sym)
    ctx.newpen(1)


def _cart(dec_deg: float, inc_deg: float) -> np.ndarray:
    """Equivalent de `cart(1.,d,ai,...)` (anisotropie.f:978-984)."""
    d = math.radians(dec_deg)
    i = math.radians(inc_deg)
    return np.array([math.cos(d) * math.cos(i), math.sin(d) * math.cos(i), math.sin(i)])


def draw_stereo_confidence_ellipses(
    ctx: PlotContext,
    results: List[TensorialMeanResult],
    r: float,
    invert_negative: bool = True,
    ams_iproj: int = 0,
) -> None:
    """Equivalent du trace des ellipses de confiance a 95% de `tratm`
    (anisotropie.f:2438-2477) : pour chaque axe k1/k2/k3 d'un tenseur
    moyen, reconstruit le repere tangent (axes j,k voisins cycliques +
    l'axe i lui-meme, via `cart` sur leurs dec/inc DEJA stockes - meme
    procede que le Fortran, pas besoin de reconserver les vecteurs propres
    d'origine), parametrise le contour de l'ellipse (`restm[i]` applique a
    73 points du cercle unite, meme pas de 5 deg), projette chaque point
    (`stecor`) et relie (`plott`).

    Croisement d'hemisphere - CORRIGE (demande explicite utilisateur :
    "dans les stereo, les ellipses sortent des stereo. Normalement, la
    partie des ellipses qui sortent de la demi-sphere inferieure sont
    projetes en pointille dans l'hemisphere superieure") : un point de
    l'ellipse dont `world[2]` (inclinaison) est negatif est projete avec
    `stecor(dec, abs(inc), ...)` - MEME azimut, rayon "replie" comme s'il
    etait dans l'hemisphere inferieur (meme principe que `tratm`, qui
    calcule `ai=atan(sqrt(z*z/(1-z*z)))` = asin(|z|) - TOUJOURS positif -
    anisotropie.f:2455) au lieu du calcul signe direct utilise
    precedemment ici, qui projetait ces points HORS du cercle (rayon >
    r). Le segment est trace en pointille (`linestyle="dashed"`) des que
    son point d'ARRIVEE est dans l'hemisphere superieur - le Fortran
    d'origine se contente d'un simple lever de crayon (`ns`) a la toute
    premiere transition puis redessine plein (pas de vrai pointille sur
    un traceur a plume) ; ici, avec un vrai style de trait disponible, le
    segment entier dans l'hemisphere superieur reste pointille plutot que
    de redevenir plein apres une seule coupure, pour rendre la portion
    "remontee" immediatement identifiable sur toute sa longueur."""
    for res in results:
        if not res.axes or any(ax.restm is None for ax in res.axes):
            continue
        for i in range(3):
            j, k = (i + 1) % 3, (i + 2) % 3
            axis_i, axis_j, axis_k = res.axes[i], res.axes[j], res.axes[k]
            tm_i, tm_j, tm_k = _cart(axis_i.dec, axis_i.inc), _cart(axis_j.dec, axis_j.inc), _cart(axis_k.dec, axis_k.inc)
            flip = invert_negative and axis_i.inc < 0.0

            ctx.newpen(_AXIS_PEN[i + 1])
            ctx.thickn((i + 1) * 0.3)
            first = True
            for it in range(1, 74):
                te = math.radians(it * 5.0)
                xi = np.array([math.cos(te), math.sin(te)])
                xf = axis_i.restm @ xi
                world = xf[0] * tm_j + xf[1] * tm_k + 1.0 * tm_i
                norm = float(np.linalg.norm(world))
                if norm > 0:
                    world = world / norm
                if flip:
                    world = -world
                dec = math.degrees(math.atan2(world[1], world[0]))
                inc = math.degrees(math.asin(max(-1.0, min(1.0, world[2]))))
                x, y = stecor(dec, abs(inc), r, ams_iproj)
                linestyle = "dashed" if world[2] < 0.0 else "solid"
                _plott(ctx, x, y, 3 if first else 2, linestyle=linestyle)
                first = False
    ctx.newpen(1)
    ctx.thickn(0.5)


def draw_stereo_mean_results(
    ctx: PlotContext,
    results: List[TensorialMeanResult],
    r: float,
    invert_negative: bool = True,
    ams_iproj: int = 0,
    point_size: float = 0.324,
) -> None:
    """Equivalent de `tratm` (anisotropie.f:2391-2481, trace des axes
    k1/k2/k3 des tenseurs moyens deja calcules + leurs ellipses de
    confiance a 95%, voir `draw_stereo_confidence_ellipses`). Meme
    convention symbole/couleur que `draw_stereo_axes`, mais plus gros
    (_MEAN_SIZE_RATIO) pour distinguer un resultat moyen d'une mesure
    individuelle - demande explicite utilisateur ("is it possible to
    reduce the size of the symbols for mean tensors") : le port utilisait
    un facteur 2.0x arbitraire, plus grand que le VRAI ratio Fortran
    (`stereo`, h=r/15 pour les moyennes vs h=r/25 pour les mesures
    individuelles, anisotropie.f:2091/2117 - soit 25/15 = 5/3)."""
    ctx.thickn(1.0)
    for res in results:
        if not res.axes:
            continue
        for j, ax in enumerate(res.axes, start=1):
            dec, inc = ax.dec, ax.inc
            sym = _IPOS[j]
            if inc < 0.0:
                inc = -inc
                if invert_negative:
                    dec += 180.0
                    if dec > 360.0:
                        dec -= 360.0
                else:
                    sym = _NEG[j]
            x, y = stecor(dec, inc, r, ams_iproj)
            ctx.newpen(_AXIS_PEN[j])
            _symbo1(ctx, x, y, point_size * _MEAN_SIZE_RATIO, sym)
    ctx.newpen(1)
    draw_stereo_confidence_ellipses(ctx, results, r, invert_negative, ams_iproj)


def build_stereo_figure(
    measurements: List[AMSMeasurement],
    orientation: int = 2,
    invert_negative: bool = True,
    ams_iproj: int = 0,
    mean_results: Optional[List[TensorialMeanResult]] = None,
    fig: Optional[Figure] = None,
) -> Figure:
    """Equivalent de `sterams`/`stereo` (menu Graphics > Stereographique) :
    cadre du reseau + axes k1(carre)/k2(triangle)/k3(cercle) de chaque
    mesure de `measurements`, en couleur (rouge/vert/bleu). `mean_results`
    (optionnel) : resultats de moyenne tensorielle deja calcules, traces en
    plus (equivalent partiel de `tratm`, symboles 2x plus grands pour les
    distinguer des mesures individuelles)."""
    dimster = 15.0  # rdim toujours 5.0*3 cote Fortran, quel que soit le prompt (bug reel confirme)
    point_size = (0.18 * dimster) / 10.0

    if fig is None:
        fig = Figure(figsize=(5.5, 5.5), dpi=100)
    else:
        fig.clear()
    ax = fig.add_subplot(111)
    ctx = PlotContext(ax)
    ctx.clear()
    ctx.plot(0.0, 0.0, -3)

    r = draw_stereo_net(ctx, orientation, ams_iproj, dimster)
    draw_stereo_axes(ctx, measurements, orientation, r, invert_negative, ams_iproj, point_size)
    if mean_results:
        draw_stereo_mean_results(ctx, mean_results, r, invert_negative, ams_iproj, point_size)

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="s", color=_AXIS_COLOR[1], linestyle="", label="k1", markerfacecolor=_AXIS_COLOR[1]),
        Line2D([0], [0], marker="^", color=_AXIS_COLOR[2], linestyle="", label="k2", markerfacecolor=_AXIS_COLOR[2]),
        Line2D([0], [0], marker="o", color=_AXIS_COLOR[3], linestyle="", label="k3", markerfacecolor=_AXIS_COLOR[3]),
    ]
    if mean_results:
        handles.append(Line2D([0], [0], marker="s", color="black", linestyle="",
                               markerfacecolor="none", markersize=8, label="mean (5/3x)"))
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=False)

    ax.relim()
    ax.autoscale_view()
    ax.set_title(", ".join(m.id for m in measurements) if len(measurements) <= 3 else f"{len(measurements)} mesures")
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------
# Bootstrap (ams_bootstrap.py) : PAS un port Fortran (bootams est un stub
# vide, voir ams_bootstrap.py) - reproduit ce que PmagPy propose deja
# (ipmag.plot_aniso, options `ivec`/`iboot`) : le nuage des directions
# rechantillonnees ("ivec=True"), et/ou l'ellipse de Kent (zeta/eta,
# "ivec=False") - meme moteur graphique (`stecor`/`_plott`/PlotContext)
# que le reste du stereonet AMS, pour rester coherent visuellement.
# ----------------------------------------------------------------------

def draw_stereo_bootstrap_cloud(
    ctx: PlotContext,
    result: BootstrapResult,
    r: float,
    invert_negative: bool = True,
    ams_iproj: int = 0,
    point_size: float = 0.08,
) -> None:
    """Equivalent de `ipmag.plot_aniso(..., iboot=True, ivec=True)` : trace
    la direction (dec,inc) de CHAQUE tirage bootstrap (petit point, pas de
    symbole rempli) pour les 3 axes - visualise directement la dispersion,
    sans passer par l'ellipse de Kent resumee."""
    ctx.thickn(0.3)
    for i, axis in enumerate(result.axes, start=1):
        ctx.newpen(_AXIS_PEN[i])
        for dec, inc in axis.draws:
            if inc < 0.0:
                inc = -inc
                if invert_negative:
                    dec += 180.0
                    if dec > 360.0:
                        dec -= 360.0
            x, y = stecor(dec, inc, r, ams_iproj)
            _symbo1(ctx, x, y, point_size, 2)  # type 2 -> remappe en triangle plein (petit)
    ctx.newpen(1)


def draw_stereo_bootstrap_ellipse(
    ctx: PlotContext,
    result: BootstrapResult,
    r: float,
    invert_negative: bool = True,
    ams_iproj: int = 0,
) -> None:
    """Equivalent de `ipmag.plot_aniso(..., iboot=True, ivec=False)` : trace
    l'ellipse de confiance de Kent (zeta/eta, `ams_bootstrap.sbootpars`)
    autour de la direction moyenne de chaque axe - meme principe
    geometrique que `draw_stereo_confidence_ellipses` (Jelinek), mais SANS
    matrice `restm` : zeta/eta sont deja les 2 demi-axes le long de leurs
    propres directions tangentes (zeta_dec/inc, eta_dec/inc), donc le
    contour est directement `(zeta*cos(theta), eta*sin(theta))` (diagonal,
    pas de terme croise)."""
    for axis in result.axes:
        tm_mean = _cart(axis.dec, axis.inc)
        tm_zeta = _cart(axis.zeta_dec, axis.zeta_inc)
        tm_eta = _cart(axis.eta_dec, axis.eta_inc)
        zeta_rad = math.radians(axis.zeta)
        eta_rad = math.radians(axis.eta)
        flip = invert_negative and axis.inc < 0.0

        ctx.newpen(1)
        ctx.thickn(0.8)
        first = True
        for it in range(1, 74):
            te = math.radians(it * 5.0)
            xf0 = zeta_rad * math.cos(te)
            xf1 = eta_rad * math.sin(te)
            world = xf0 * tm_zeta + xf1 * tm_eta + 1.0 * tm_mean
            norm = float(np.linalg.norm(world))
            if norm > 0:
                world = world / norm
            if flip:
                world = -world
            dec = math.degrees(math.atan2(world[1], world[0]))
            inc = math.degrees(math.asin(max(-1.0, min(1.0, world[2]))))
            x, y = stecor(dec, inc, r, ams_iproj)
            _plott(ctx, x, y, 3 if first else 2)
            first = False
    ctx.newpen(1)
    ctx.thickn(0.5)


def build_bootstrap_stereo_figure(
    result: BootstrapResult,
    orientation: int = 2,
    invert_negative: bool = True,
    ams_iproj: int = 0,
    show_cloud: bool = True,
    show_ellipse: bool = True,
    fig: Optional[Figure] = None,
) -> Figure:
    """Equivalent (reseau AMS + options PmagPy `ivec`/`iboot`) du 2e
    graphique de `ipmag.plot_aniso` (moyenne + confiance bootstrap) : nuage
    des tirages et/ou ellipse de Kent, superposes a la direction moyenne."""
    dimster = 15.0
    point_size = (0.18 * dimster) / 10.0

    if fig is None:
        fig = Figure(figsize=(5.5, 5.5), dpi=100)
    else:
        fig.clear()
    ax = fig.add_subplot(111)
    ctx = PlotContext(ax)
    ctx.clear()
    ctx.plot(0.0, 0.0, -3)

    r = draw_stereo_net(ctx, orientation, ams_iproj, dimster, show_orient_label=False)
    if show_cloud:
        draw_stereo_bootstrap_cloud(ctx, result, r, invert_negative, ams_iproj)
    if show_ellipse:
        draw_stereo_bootstrap_ellipse(ctx, result, r, invert_negative, ams_iproj)
    # direction moyenne, symboles pleins standard (carre/triangle/cercle)
    for i, axis in enumerate(result.axes, start=1):
        dec, inc = axis.dec, axis.inc
        sym = _IPOS[i]
        if inc < 0.0:
            inc = -inc
            if invert_negative:
                dec += 180.0
                if dec > 360.0:
                    dec -= 360.0
            else:
                sym = _NEG[i]
        x, y = stecor(dec, inc, r, ams_iproj)
        ctx.newpen(_AXIS_PEN[i])
        _symbo1(ctx, x, y, point_size * 1.5, sym)
    ctx.newpen(1)

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="s", color=_AXIS_COLOR[1], linestyle="", label="k1", markerfacecolor=_AXIS_COLOR[1]),
        Line2D([0], [0], marker="^", color=_AXIS_COLOR[2], linestyle="", label="k2", markerfacecolor=_AXIS_COLOR[2]),
        Line2D([0], [0], marker="o", color=_AXIS_COLOR[3], linestyle="", label="k3", markerfacecolor=_AXIS_COLOR[3]),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=False)

    ax.relim()
    ax.autoscale_view()
    method = "parametrique" if result.parametric else "non parametrique"
    ax.set_title(f"Bootstrap ({method}, nb={result.nb}, n={result.n})")
    fig.tight_layout()
    return fig
