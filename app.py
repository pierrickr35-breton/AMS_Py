"""
AMS_Py : port Python de AMS_OSX_AWE (reference/AMS_OSX_AWE), application
Fortran/Qt5 de traitement d'anisotropie de susceptibilite magnetique (AMS),
meme heritage de code et meme logique d'interface que STARpaleomag_Py (console
texte + panneau graphique, menu Fortran reproduit tel quel).

ETAT (premiere passe) : fichiers .ANI (ouverture/liste), selection/liste des
mesures, orientation (echantillon/in-situ/pendage corrige), et deux routines
de Calcul entierement portees et verifiees (Inverse correction, Substraction
of tensors). Tensorial mean/Mean Results/Graphics dependent des statistiques
de Hext/Jelinek (anisotropie.f, tsmean/fisherams/flinn/tpprim/stereo) - EN
COURS DE PORTAGE, stubs clairement marques ci-dessous."""

import math
import os
import subprocess
import sys
import webbrowser
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ams_selection import (
    AMSMeasurement,
    read_ani_file,
    select_measurements,
    list_measurements,
    apply_orientation,
    import_legacy_ani,
    mean_result_site,
    write_ani_mean_result,
    read_ani_mean_results,
    read_ani_mean_results_from_pmagani,
    is_imaginary_component,
    mark_pmagani_export,
    create_empty_pmagani_if_missing,
    _ORIENT_TO_FILE_CODE,
)
from ams_prmag import read_prmag_specimens, ani_path_for, create_prmag_from_legacy_ani, create_prmag_from_asc
from ams_asc import archive_asc_file
from ams_calcul import correct_direction_with_tensor, subtract_tensors
from ams_bootstrap import compute_bootstrap_mean, format_bootstrap_result
from ams_stats import (
    tsmean, mean_susceptibility, format_measurement_list, format_tsmean_box,
    format_lisresmem_table, magic_site_aniso_fields, principal_axes,
    negative_k_corrected_tensor, shape_params,
)
from ams_stereo import build_stereo_figure, build_bootstrap_stereo_figure
from ams_xy import (
    build_anisotropy_parameters_figure,
    shape_entries_from_measurements,
    shape_entries_from_mean_results,
    build_susceptibility_anisotropy_figure,
    build_im_re_susceptibility_figure,
    P_OUTLIER_THRESHOLD,
)


# Raccourcis clavier : les libelles de menu affichaient deja "(Cmd+X)" mais
# aucun binding reel n'existait ("les commandes ne sont pas activees dans
# les menus") - meme convention que STARpaleomag_Py/app.py (Command comme
# modificateur de base, Cmd+Ctrl pour les combos "META"/"ALT" du Fortran,
# jamais Option/Alt seul - piege AZERTY documente cote STARpaleomag_Py). Chaque
# entree est (texte affiche, sequence Tk) - voir _labeled.
_SHORTCUTS_MAC = {
    "openani": ("Cmd+O",      "<Command-o>"),
    "selmes":  ("Cmd+E",      "<Command-e>"),
    "lismes":  ("Cmd+L",      "<Command-l>"),
    "initmes": ("Cmd+I",      "<Command-i>"),
    "delmes":  ("Cmd+D",      "<Command-d>"),
    "tsmean":  ("Cmd+M",      "<Command-m>"),
    "stereo":  ("Cmd+Ctrl+S", "<Command-Control-s>"),
    "aniso":   ("Cmd+Ctrl+P", "<Command-Control-p>"),
    # Coordonnees (sample/in-situ/tilt corrected) - MEMES combinaisons que
    # STARpaleomag_Py/app.ORIENTATION_SHORTCUT_NAMES ("selce"/"selis"/"selcp",
    # Cmd+Ctrl+A/B/T) - demande explicite utilisateur ("can you provide in
    # AMS_Py the same shortcuts as in STARpaleomag_Py for coordinates"). "aniso"
    # (Anisotropy parameters) a du etre deplace de Cmd+Ctrl+A vers
    # Cmd+Ctrl+P pour liberer le "A" au profit de "selce", comme cote
    # STARpaleomag_Py.
    "selce":   ("Cmd+Ctrl+A", "<Command-Control-a>"),
    "selis":   ("Cmd+Ctrl+B", "<Command-Control-b>"),
    "selcp":   ("Cmd+Ctrl+T", "<Command-Control-t>"),
}

# Version PC/Windows - meme raisonnement que STARpaleomag_Py/app._SHORTCUTS_WIN
# (demande explicite utilisateur "provide shortcuts for both mac and
# windows") : "Command"/"Meta" ne correspond a aucune touche physique
# fiable sous Tk sur Windows, et Ctrl+lettre seul percute les bindings
# Emacs internes de Tk (Entry/Text) ainsi que les raccourcis Windows
# standard - Ctrl+Shift+<LETTRE> est donc le seul palier retenu, PAS de
# distinction Cmd+lettre/Cmd+Ctrl+lettre comme sur Mac (chaque commande
# garde neanmoins la MEME lettre qu'en Mac : les 11 lettres du dict Mac
# sont deja toutes distinctes d'un palier a l'autre, donc aucune
# reaffectation n'est necessaire en les regroupant ici).
_SHORTCUTS_WIN = {
    "openani": ("Ctrl+Shift+O", "<Control-Shift-O>"),
    "selmes":  ("Ctrl+Shift+E", "<Control-Shift-E>"),
    "lismes":  ("Ctrl+Shift+L", "<Control-Shift-L>"),
    "initmes": ("Ctrl+Shift+I", "<Control-Shift-I>"),
    "delmes":  ("Ctrl+Shift+D", "<Control-Shift-D>"),
    "tsmean":  ("Ctrl+Shift+M", "<Control-Shift-M>"),
    "stereo":  ("Ctrl+Shift+S", "<Control-Shift-S>"),
    "aniso":   ("Ctrl+Shift+P", "<Control-Shift-P>"),
    "selce":   ("Ctrl+Shift+A", "<Control-Shift-A>"),
    "selis":   ("Ctrl+Shift+B", "<Control-Shift-B>"),
    "selcp":   ("Ctrl+Shift+T", "<Control-Shift-T>"),
}

SHORTCUTS = _SHORTCUTS_WIN if sys.platform.startswith("win") else _SHORTCUTS_MAC


def _resource_path(*parts: str) -> str:
    """Chemin absolu d'une ressource livree AVEC l'appli (le guide
    utilisateur HTML, voir ouvrir_user_guide) - equivalent de
    STARpaleomag_Py/app._resource_path : fonctionne aussi bien lancee
    depuis les sources (repertoire de ce fichier) qu'empaquetee par
    PyInstaller (sys._MEIPASS - voir AMS_Py.spec, datas += [('help',
    'help')]) - demande explicite utilisateur ("is it possible to write a
    guide for AMS_Py")."""
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, *parts)


class AmsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("AMS_Py - Anisotropy of Magnetic Susceptibility")
        self.root.geometry("1100x700")

        self.donnees: list = []    # toutes les mesures du fichier .ANI (equiv amsfichier)
        self.selection: list = []  # mesures selectionnees (equiv ams(:))
        self.ani_path = None
        self.prmag_path = None
        self.prmag_specimens = {}  # {specimen_id: PrmagSpecimen} - site/sample MagIC, voir ams_prmag
        self.orientation = tk.IntVar(value=2)  # 1=echantillon,2=in-situ,3=pendage corrige (defaut in-situ, comme iorient=2)
        # Ancre du guide utilisateur correspondant au DERNIER item de menu
        # invoque (voir _menu_cmd/ouvrir_user_guide, meme mecanisme que
        # STARpaleomag_Py/app.py) - None tant qu'aucun menu de contenu n'a
        # encore ete utilise, auquel cas le guide s'ouvre sur sa page
        # d'accueil.
        self._help_anchor = None
        self.mean_results = []  # [(orientation, TensorialMeanResult), ...] - equivalent minimal de amsres (Mean Results, pas encore porte en menu complet)
        self._current_graphic = None
        self.ams_iproj = 0  # 'paramster' : 0=lambert/equiaire (defaut Fortran), 1=stereographique
        self.invert_negative = True  # 'iinv' (paramster) : defaut 'y'
        self._aniso_source = 0  # 'Anisotropy parameters' : 0=donnees,1=resultats,2=d+r
        # bornes hautes d'echelle (haut/bas) choisies par l'utilisateur
        # pour "Anisotropy parameters" - None = auto (voir ams_xy.build_
        # anisotropy_parameters_figure) - demande explicite utilisateur
        # ("est-ce possible de selectionner l'echelle"), INDEPENDANTES
        # entre les deux panneaux (correction ulterieure : "the scale
        # should not apply to both plots (Re et Im) as they have very
        # different ranges").
        self._aniso_scale_top = None
        self._aniso_scale_bottom = None
        self._last_bootstrap = None  # dernier BootstrapResult calcule (pour le plot bootstrap)
        self._bootstrap_show_cloud = True
        self._bootstrap_show_ellipse = True

        self._setup_menu()
        self._setup_shortcuts()

        self.paned_window = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        self.paned_window.pack(fill=tk.BOTH, expand=True)

        self.graph_frame = ttk.Frame(self.paned_window, width=550)
        self.fig = Figure(figsize=(5.2, 5.2), dpi=100)
        self.canvas_fig = FigureCanvasTkAgg(self.fig, master=self.graph_frame)
        self.canvas_fig.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.paned_window.add(self.graph_frame, weight=1)
        # Clic-pour-info sur un point du graphique "Anisotropy parameters"
        # (id specimen) et clic-sur-axe-X pour en changer l'echelle - voir
        # ams_xy._register_point_pick/_register_xaxis_pick - demande
        # explicite utilisateur ("click on point to have the specimen
        # number" / "click on the X axis to change the max value").
        self.canvas_fig.mpl_connect("pick_event", self._on_plot_pick)

        self.text_frame = ttk.Frame(self.paned_window, width=550)
        self.text_area = tk.Text(
            self.text_frame, bg="#ffffff", fg="#000000", insertbackground="black",
            font=("Courier", 14), wrap="none",
        )
        text_yscroll = ttk.Scrollbar(self.text_frame, orient=tk.VERTICAL, command=self.text_area.yview)
        text_xscroll = ttk.Scrollbar(self.text_frame, orient=tk.HORIZONTAL, command=self.text_area.xview)
        self.text_area.configure(yscrollcommand=text_yscroll.set, xscrollcommand=text_xscroll.set)
        self.text_area.tag_configure("prompt", foreground="#c0392b")
        # Verdict PmagPy Hext F-test (satisfactory/not satisfactory) mis
        # en evidence - demande explicite utilisateur ("write the message
        # in bold green when it is satisfactory and in bold red when it
        # is not, especially in the view the 15 ATRM tensors per
        # sample") - voir _afficher.
        self.text_area.tag_configure("sat_ok", foreground="#1e8449", font=("Courier", 14, "bold"))
        self.text_area.tag_configure("sat_bad", foreground="#c0392b", font=("Courier", 14, "bold"))
        self.text_area.grid(row=0, column=0, sticky="nsew")
        text_yscroll.grid(row=0, column=1, sticky="ns")
        text_xscroll.grid(row=1, column=0, sticky="ew")
        self.text_frame.rowconfigure(0, weight=1)
        self.text_frame.columnconfigure(0, weight=1)
        self.paned_window.add(self.text_frame, weight=1)

        self.root.after(200, self._activate_window)

    def _activate_window(self):
        self.root.attributes("-topmost", True)
        self.root.after(50, lambda: self.root.attributes("-topmost", False))
        try:
            subprocess.run(
                ["osascript", "-e",
                 f'tell application "System Events" to set frontmost of '
                 f'(first process whose unix id is {os.getpid()}) to true'],
                check=False, capture_output=True, timeout=2,
            )
        except (subprocess.SubprocessError, OSError):
            pass

    # ------------------------------------------------------------------
    # Console texte (meme motif que STARpaleomag_Py/app.py)
    # ------------------------------------------------------------------

    def _afficher(self, text):
        if self.text_area.get("1.0", "end-1c").strip():
            self.text_area.insert(tk.END, "\n" + "-" * 60 + "\n")
        # Met en gras vert/rouge toute ligne portant le verdict PmagPy
        # Hext F-test (satisfactory/not satisfactory) - demande explicite
        # utilisateur ("write the message in bold green when it is
        # satisfactory and in bold red when it is not, especially in the
        # view the 15 ATRM tensors per sample") - voir ams_stats.
        # format_measurement_list, qui produit cette ligne pour chaque
        # groupe de specimen affiche par "List measurements"/etc.
        for line in text.splitlines(keepends=True):
            if "(not satisfactory)" in line:
                self.text_area.insert(tk.END, line, "sat_bad")
            elif "(satisfactory)" in line:
                self.text_area.insert(tk.END, line, "sat_ok")
            else:
                self.text_area.insert(tk.END, line)
        self.text_area.see(tk.END)

    def _console_input(self, prompt, default=""):
        self.text_area.insert(tk.END, prompt, "prompt")
        start_index = self.text_area.index("end-1c")
        self.text_area.mark_set(tk.INSERT, tk.END)
        self.text_area.see(tk.END)
        self.text_area.focus_set()

        outcome = {"value": None}
        done = tk.BooleanVar(value=False)

        def on_return(event):
            typed = self.text_area.get(start_index, "end-1c")
            outcome["value"] = typed if typed.strip() else default
            self.text_area.insert(tk.END, "\n")
            done.set(True)
            return "break"

        def on_escape(event):
            self.text_area.delete(start_index, "end-1c")
            outcome["value"] = None
            self.text_area.insert(tk.END, "\n")
            done.set(True)
            return "break"

        ret_id = self.text_area.bind("<Return>", on_return)
        esc_id = self.text_area.bind("<Escape>", on_escape)
        self.text_area.wait_variable(done)
        self.text_area.unbind("<Return>", ret_id)
        self.text_area.unbind("<Escape>", esc_id)
        self.text_area.see(tk.END)
        return outcome["value"]

    def _restore_focus(self):
        """Apres un popup natif (messagebox), le focus clavier ne revient
        PAS automatiquement sur la fenetre principale (constate sur macOS/
        Aqua) - les raccourcis lies via bind_all restent inertes tant qu'on
        n'a pas clique manuellement dans une fenetre de l'appli - demande
        explicite utilisateur ("we need to manually click in one of the
        window before to activate the shortcuts. Is there a way to change
        this behavior?"). Force le focus sur la fenetre principale, et sur
        la console texte (cible naturelle de la frappe/des raccourcis)."""
        self.root.lift()
        self.root.focus_force()
        self.text_area.focus_set()

    def _showinfo(self, *args, **kwargs):
        result = messagebox.showinfo(*args, **kwargs)
        self._restore_focus()
        return result

    def _showwarning(self, *args, **kwargs):
        result = messagebox.showwarning(*args, **kwargs)
        self._restore_focus()
        return result

    def _showerror(self, *args, **kwargs):
        result = messagebox.showerror(*args, **kwargs)
        self._restore_focus()
        return result

    def _not_implemented(self, name):
        self._showinfo(
            "Not yet ported",
            f"'{name}' is not yet ported to Python (pending anisotropie.f "
            f"statistics - Hext/Jelinek tensorial mean).",
        )

    # ------------------------------------------------------------------
    # Menu
    # ------------------------------------------------------------------

    @staticmethod
    def _labeled(text, shortcut_name):
        """Libelle avec le raccourci entre parentheses, SANS `accelerator=`
        (voir STARpaleomag_Py/app._labeled - meme raison : sur Aqua, `accelerator=`
        semble faire intercepter la combinaison par le menu natif sans
        jamais invoquer la commande Tcl, court-circuitant bind_all)."""
        return f"{text}    ({SHORTCUTS[shortcut_name][0]})"

    def _menu_cmd(self, anchor, func):
        """Enveloppe `func` (la vraie commande d'un item de menu) pour
        enregistrer `anchor` comme dernier contexte d'aide AVANT de
        l'executer - meme mecanisme que STARpaleomag_Py/app._menu_cmd.
        `anchor` correspond a l'id du BLOC (separateur a separateur) dans
        help/AMS_Py_Guide.html ; plusieurs items d'un meme bloc partagent
        la meme ancre."""
        def wrapped(*args, **kwargs):
            self._help_anchor = anchor
            return func(*args, **kwargs)
        return wrapped

    def _setup_menu(self):
        menubar = tk.Menu(self.root)

        files_menu = tk.Menu(menubar, tearoff=0)
        files_menu.add_command(
            label=self._labeled("Open File .pmagani...", "openani"),
            command=self._menu_cmd("files-open", self.ouvrir_ani_dialog))
        files_menu.add_command(
            label="Open File .prmag...", command=self._menu_cmd("files-open", self.ouvrir_prmag_dialog))
        files_menu.add_command(
            label="List File .pmagani", command=self._menu_cmd("files-open", self.lister_fichier_ani))
        files_menu.add_separator()
        files_menu.add_command(
            label="Import legacy .ANI to .pmagani...",
            command=self._menu_cmd("files-legacy", self.ouvrir_import_legacy_ani_dialog))
        files_menu.add_command(
            label="Create prmag & pmagani from legacy .ANI...",
            command=self._menu_cmd("files-legacy", self.ouvrir_creer_prmag_from_ani_dialog))
        files_menu.add_command(
            label="Create prmag & pmagani from .asc AGICO file...",
            command=self._menu_cmd("files-legacy", self.ouvrir_creer_prmag_from_asc_dialog))
        files_menu.add_separator()
        files_menu.add_command(
            label="Archive ASC into .pmagani...",
            command=self._menu_cmd("files-archive", self.ouvrir_archiver_asc_dialog))
        files_menu.add_separator()
        files_menu.add_command(
            label="Mark selection for MagIC export...",
            command=self._menu_cmd("files-export", self.ouvrir_marquer_export_dialog))
        menubar.add_cascade(label="AMS Files", menu=files_menu)

        data_menu = tk.Menu(menubar, tearoff=0)
        data_menu.add_command(
            label=self._labeled("Select measurements", "selmes"),
            command=self._menu_cmd("data-list", self.selectionner_mesures))
        data_menu.add_command(
            label=self._labeled("List measurements", "lismes"),
            command=self._menu_cmd("data-list", self.lister_mesures))
        data_menu.add_command(
            label="List measurements with depth",
            command=self._menu_cmd("data-list", lambda: self._not_implemented("List measurements with depth")))
        data_menu.add_separator()
        data_menu.add_command(
            label=self._labeled("init list to zero", "initmes"),
            command=self._menu_cmd("data-edit", self.reinitialiser_liste))
        data_menu.add_command(
            label=self._labeled("Delete lines", "delmes"),
            command=self._menu_cmd("data-edit", self.supprimer_lignes))
        data_menu.add_separator()
        data_menu.add_radiobutton(
            label=self._labeled("Sample coordinates", "selce"), variable=self.orientation, value=1,
            command=self._menu_cmd("data-orient", self.changer_orientation))
        data_menu.add_radiobutton(
            label=self._labeled("In situ coordinates", "selis"), variable=self.orientation, value=2,
            command=self._menu_cmd("data-orient", self.changer_orientation))
        data_menu.add_radiobutton(
            label=self._labeled("Tilt corrected coordinates", "selcp"), variable=self.orientation, value=3,
            command=self._menu_cmd("data-orient", self.changer_orientation))
        menubar.add_cascade(label="AMS data", menu=data_menu)

        results_menu = tk.Menu(menubar, tearoff=0)
        results_menu.add_command(
            label="Select results", command=self._menu_cmd("results-manage", self.ouvrir_selres_dialog))
        results_menu.add_command(
            label="List results", command=self._menu_cmd("results-manage", self.ouvrir_lisresmem_dialog))
        results_menu.add_command(
            label="init list results", command=self._menu_cmd("results-manage", self.ouvrir_initres_dialog))
        results_menu.add_command(
            label="delete result", command=self._menu_cmd("results-manage", self.ouvrir_delres_dialog))
        results_menu.add_command(
            label="Liste results file", command=self._menu_cmd("results-manage", self.ouvrir_lisresfich_dialog))
        results_menu.add_command(
            label="save results in File", command=self._menu_cmd("results-manage", self.ouvrir_sauveresfich_dialog))
        menubar.add_cascade(label="Mean Results", menu=results_menu)

        calcul_menu = tk.Menu(menubar, tearoff=0)
        calcul_menu.add_command(
            label=self._labeled("Tensorial mean", "tsmean"),
            command=self._menu_cmd("calcul-tsmean", self.ouvrir_tsmean_dialog))
        calcul_menu.add_separator()
        calcul_menu.add_command(
            label="Mean susceptibility", command=self._menu_cmd("calcul-msus", self.ouvrir_mds_dialog))
        calcul_menu.add_separator()
        calcul_menu.add_command(
            label="Fisher 3 axes",
            command=self._menu_cmd("calcul-fisher", lambda: self._not_implemented("Fisher 3 axes")))
        calcul_menu.add_command(
            label="ellipses Bootstrap", command=self._menu_cmd("calcul-fisher", self.ouvrir_bootstrap_dialog))
        calcul_menu.add_separator()
        calcul_menu.add_command(
            label="Inverse correction", command=self._menu_cmd("calcul-corr", self.ouvrir_corpaldir_dialog))
        calcul_menu.add_command(
            label="Substraction of tensors", command=self._menu_cmd("calcul-corr", self.ouvrir_soustract_dialog))
        menubar.add_cascade(label="Calcul", menu=calcul_menu)

        graphics_menu = tk.Menu(menubar, tearoff=0)
        graphics_menu.add_command(
            label=self._labeled("Stereographique", "stereo"),
            command=self._menu_cmd("graphics-stereo", self.afficher_stereo))
        graphics_menu.add_command(
            label="Stereographic sample",
            command=self._menu_cmd("graphics-stereo", self.ouvrir_stereo_sample_dialog))
        graphics_menu.add_command(
            label="Stereographic site", command=self._menu_cmd("graphics-stereo", self.ouvrir_stereo_site_dialog))
        graphics_menu.add_command(
            label="parametres projection", command=self._menu_cmd("graphics-stereo", self.ouvrir_paramster_dialog))
        graphics_menu.add_separator()
        graphics_menu.add_command(
            label=self._labeled("Anisotropy parameters", "aniso"),
            command=self._menu_cmd("graphics-aniso", self.afficher_anisotropy_parameters))
        graphics_menu.add_separator()
        graphics_menu.add_command(
            label="Susceptibility / anisotropy degree",
            command=self._menu_cmd("graphics-shape", self.afficher_susceptibility_anisotropy))
        graphics_menu.add_command(
            label="Im vs Re susceptibility",
            command=self._menu_cmd("graphics-shape", self.afficher_im_re_susceptibility))
        graphics_menu.add_separator()
        graphics_menu.add_command(
            label="Export SVG...", command=self._menu_cmd("graphics-svg", self.exporter_svg))
        menubar.add_cascade(label="Graphics", menu=graphics_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="User Guide", command=self.ouvrir_user_guide)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.root.config(menu=menubar)

    def _set_orientation(self, value):
        """Equivalent de STARpaleomag_Py/app._set_orientation : bascule
        self.orientation puis rejoue changer_orientation - utilise par les
        raccourcis clavier selce/selis/selcp (les radiobuttons du menu le
        font deja via leur propre `command`, mais un raccourci clavier
        n'active pas le radiobutton lui-meme, juste sa `command`)."""
        self.orientation.set(value)
        self.changer_orientation()

    def _setup_shortcuts(self):
        """Lie reellement les raccourcis affiches dans les libelles de menu
        (voir SHORTCUTS) - jusqu'ici purement decoratifs."""
        bindings = {
            "openani": self.ouvrir_ani_dialog,
            "selmes": self.selectionner_mesures,
            "lismes": self.lister_mesures,
            "initmes": self.reinitialiser_liste,
            "delmes": self.supprimer_lignes,
            "tsmean": self.ouvrir_tsmean_dialog,
            "stereo": self.afficher_stereo,
            "aniso": self.afficher_anisotropy_parameters,
            "selce": lambda: self._set_orientation(1),
            "selis": lambda: self._set_orientation(2),
            "selcp": lambda: self._set_orientation(3),
        }
        for name, callback in bindings.items():
            self.root.bind_all(SHORTCUTS[name][1], lambda event, cb=callback: cb())

    # ------------------------------------------------------------------
    # AMS Files
    # ------------------------------------------------------------------

    @staticmethod
    def _prmag_path_for(pmagani_path):
        """Inverse de ams_prmag.ani_path_for : le .prmag compagnon d'un
        .pmagani, MEME nom de base ("the user should give the same name to
        the pmagani than the .prmag file" - demande explicite
        utilisateur)."""
        return os.path.splitext(pmagani_path)[0] + ".prmag"

    def _join_prmag_orientation(self):
        """Complete cin/caz/dip/str_ de self.donnees depuis self.
        prmag_specimens (jointure par numero de specimen) - cin/caz/dip/str
        ne sont plus stockes dans .pmagani lui-meme ("the codes cin caz dip
        str are now not needed if the pmagani is linked to the prmag
        file"). Retourne la liste des id de MESURE (un par mesure sans
        correspondance dans le .prmag, PAS dedupliquee - un id peut
        apparaitre deux fois si Re et Im sont tous deux sans specimen -
        `len(...)` reste donc le meme compte "N measurement(s)" qu'avant)
        dont l'orientation reste a 0.0 - demande explicite utilisateur
        ("can you list the specimens without information") : chaque
        appelant n'affichait jusqu'ici qu'un compte, pas les id concernes
        (voir _format_missing_prmag_warning pour l'affichage, deduplique
        celui-la)."""
        missing_ids = []
        for m in self.donnees:
            spec = self.prmag_specimens.get(m.id)
            if spec is None:
                missing_ids.append(m.id)
                continue
            m.cin, m.caz, m.dip, m.str_ = spec.cin, spec.caz, spec.dip, spec.str_
        return missing_ids

    def _format_missing_prmag_warning(self, missing_ids):
        """Formate l'avertissement "N measurement(s) have no matching
        specimen in the .prmag file", en listant les id CONCERNES
        (deduplique pour l'affichage - le compte N, lui, reste celui de
        `missing_ids` tel quel, une entree par MESURE) - plafonne a 20 id
        affiches ("..." au-dela) - factorise entre les 3 appelants de
        _join_prmag_orientation."""
        unique_ids = list(dict.fromkeys(missing_ids))
        shown = ", ".join(unique_ids[:20])
        more = ", ..." if len(unique_ids) > 20 else ""
        return (
            f"warning: {len(missing_ids)} measurement(s) have no matching specimen in the "
            f".prmag file - their orientation (cin/caz/dip/str) stays at 0.0: {shown}{more}\n")

    def ouvrir_ani_dialog(self):
        """Equivalent de `openani` (lect_asc.f:592-676) - ouvre un fichier
        .pmagani (nouveau format) ou un ancien .ANI (les deux restent
        lisibles directement, voir ams_selection.read_ani_file, dispatch
        par extension). Si un .pmagani est ouvert directement ICI (plutot
        que via "Open File .prmag..."), le .prmag compagnon de MEME NOM DE
        BASE est recherche et charge automatiquement - demande explicite
        utilisateur ("the user should give the same name to the pmagani
        than the .prmag file; then opening pmagani would automatically
        open the prmag file. This will be already the case for ATRM" :
        cote STARpaleomag_Py, ouvrir le .prmag charge deja automatiquement son
        .pmagani ; ceci est la reciproque, ouvrir le .pmagani charge son
        .prmag). Sans .prmag compagnon, cin/caz/dip/str restent a 0.0 -
        avertissement affiche dans ce cas."""
        path = filedialog.askopenfilename(
            title="Open File .pmagani",
            filetypes=[("AMS .pmagani/.ANI", "*.pmagani *.ANI *.ani"), ("All files", "*.*")])
        if not path:
            return
        self.donnees = read_ani_file(path)
        self.ani_path = path
        self.selection = []
        self._afficher(f"nb measurements in file: {len(self.donnees)}\n")
        if os.path.splitext(path)[1].lower() == ".pmagani":
            prmag_candidate = self._prmag_path_for(path)
            if os.path.exists(prmag_candidate):
                self.prmag_specimens = read_prmag_specimens(prmag_candidate)
                self.prmag_path = prmag_candidate
                missing = self._join_prmag_orientation()
                self._afficher(f"companion prmag file found: {prmag_candidate}\n"
                                f"nb specimens in prmag file: {len(self.prmag_specimens)}\n")
                if missing:
                    self._afficher(self._format_missing_prmag_warning(missing))
            else:
                self._afficher(
                    "note: orientation (cin/caz/dip/str) is not stored in .pmagani anymore, and "
                    f"no companion .prmag file was found ({prmag_candidate}) - it stays at 0.0 "
                    "unless you open the matching .prmag file (AMS Files > Open File .prmag...).\n")

    def ouvrir_prmag_dialog(self):
        """Ouvre un fichier .prmag (STARpaleomag_Py) et charge automatiquement le
        .pmagani compagnon (meme nom de base - voir ams_prmag.ani_path_for,
        equivalent de calcul.ani_path_for cote STARpaleomag_Py) - demande
        explicite utilisateur ("in AMS_py we also open .prmag but it
        locates the .ANI"). Le .pmagani fournit le site/sample MagIC par
        specimen, utilises par "Export to Magic" (le .pmagani lui-meme ne
        garde que le numero de specimen - "we just keep the specimen
        number in .ani"). Si aucun .pmagani n'existe encore, retombe sur
        un ancien .ANI de meme nom de base (retro-compatibilite - voir
        "Import legacy .ANI to .pmagani..." pour le convertir)."""
        path = filedialog.askopenfilename(
            title="Open File .prmag", filetypes=[("STARpaleomag_Py .prmag", "*.prmag"), ("All files", "*.*")])
        if not path:
            return
        self.prmag_specimens = read_prmag_specimens(path)
        self.prmag_path = path
        self._afficher(f"nb specimens in prmag file: {len(self.prmag_specimens)}\n")

        candidate = ani_path_for(path)
        is_pmagani = os.path.exists(candidate)
        if not is_pmagani:
            legacy = os.path.splitext(path)[0] + ".ANI"
            if os.path.exists(legacy):
                candidate = legacy
        if os.path.exists(candidate):
            self.donnees = read_ani_file(candidate)
            self.ani_path = candidate
            self.selection = []
            missing = []
            if is_pmagani:
                # cin/caz/dip/str_ ne sont plus dans le .pmagani lui-meme
                # ("the codes cin caz dip str are now not needed if the
                # pmagani is linked to the prmag file") - on les retrouve
                # ici depuis le .prmag qu'on vient de charger, par
                # jointure sur le numero de specimen.
                missing = self._join_prmag_orientation()
            self._afficher(f"companion file found: {candidate}\n"
                            f"nb measurements in file: {len(self.donnees)}\n")
            if missing:
                self._afficher(self._format_missing_prmag_warning(missing))
        else:
            # Ni .pmagani ni ancien .ANI : en cree un vide - demande
            # explicite utilisateur ("dans AMS_Py quand on ouvre un
            # .prmag, si il n'y a pas de pmagani, en creer un vide"),
            # meme principe que STARpaleomag_Py/field_notes.
            # write_prmag_from_field_notes (compagnon cree des la
            # creation du .prmag) - ici a l'OUVERTURE cote AMS_Py, pour
            # le cas d'un .prmag deja existant (cree par une version
            # anterieure, ou par une autre voie) sans .pmagani associe.
            # `candidate` est encore ani_path_for(path) a ce stade (le
            # seul cas ou on atteint ce else) - devient directement le
            # fichier actif (self.ani_path), pret pour "Archive ASC into
            # .pmagani..." sans dialogue de fichier supplementaire.
            create_empty_pmagani_if_missing(candidate)
            self.donnees = []
            self.ani_path = candidate
            self.selection = []
            self._afficher(f"no companion .pmagani/.ANI file found - created an empty one: {candidate}\n")

    def ouvrir_import_legacy_ani_dialog(self):
        """Convertit un ANCIEN fichier .ANI (format list-directed Fortran,
        sans entete) vers le nouveau format .pmagani (tabule, avec entete
        de colonnes) - demande explicite utilisateur ("can we import old
        style .ANI in these pmagani style"), equivalent de STARpaleomag_Py/
        app.ouvrir_import_legacy_ani_dialog."""
        path = filedialog.askopenfilename(
            title="Import legacy .ANI", filetypes=[("Legacy .ANI", "*.ANI *.ani"), ("All files", "*.*")])
        if not path:
            return
        default_new_path = ani_path_for(path)
        new_path = filedialog.asksaveasfilename(
            title="Save as .pmagani", initialfile=os.path.basename(default_new_path),
            initialdir=os.path.dirname(default_new_path),
            defaultextension=".pmagani", filetypes=[("pmagani", "*.pmagani"), ("All files", "*.*")])
        if not new_path:
            return
        try:
            import_legacy_ani(path, new_path)
        except OSError as e:
            self._showerror("Import failed", f"Could not import {path}:\n{e}")
            return
        self._afficher(f"Imported legacy {path} -> {new_path}\n")
        # Compte les moyennes de site ('R' lines) importees en meme temps
        # (voir import_legacy_ani) - demande explicite utilisateur ("can
        # we check the new import legacy files where the IS TC will also
        # be written from 0 to 100") : confirme a l'utilisateur qu'elles
        # ont bien ete reportees, pas seulement les mesures brutes.
        legacy_means = read_ani_mean_results(path)
        kept = sum(1 for iorient, _r in legacy_means if iorient in (1, 2, 3))
        skipped = len(legacy_means) - kept
        if legacy_means:
            msg = f"{kept} site mean tensor result(s) also imported (with their original Sa/IS/TC orientation).\n"
            if skipped:
                msg += f"({skipped} skipped: unrecognized orientation code.)\n"
            self._afficher(msg)
        self.donnees = read_ani_file(new_path)
        self.ani_path = new_path
        self.selection = []
        prmag_candidate = self._prmag_path_for(new_path)
        if os.path.exists(prmag_candidate):
            self.prmag_specimens = read_prmag_specimens(prmag_candidate)
            self.prmag_path = prmag_candidate
            missing = self._join_prmag_orientation()
            self._afficher(f"companion prmag file found: {prmag_candidate}\n")
            if missing:
                self._afficher(self._format_missing_prmag_warning(missing))

    def ouvrir_creer_prmag_from_ani_dialog(self):
        """Cree un .prmag a partir des informations d'orientation (cin/
        caz/dip/str_) d'un ANCIEN fichier .ANI - demande explicite
        utilisateur ("un collegue souhaite avoir la possibilite de creer
        le prmag a partir du .ANI qui contient les infos de corrections
        de carotte"), equivalent AMS_Py de STARpaleomag_Py/field_notes.
        write_prmag_from_field_notes (mais depuis un .ANI plutot que des
        field notes brutes - voir ams_prmag.create_prmag_from_legacy_ani
        pour le detail, notamment la deduplication par specimen : un
        .ANI reel repete chaque specimen ATRM une fois par variante
        jackknife A0/A+/A-/A1/B1.../B6).

        Convertit AUSSI ce meme .ANI en .pmagani dans la foulee (le
        tenseur y est deja disponible, pas de raison de laisser un
        compagnon vide comme pour un .prmag cree depuis de simples
        field notes - voir ams_selection.import_legacy_ani), sauf si un
        .pmagani portant deja ce nom existe (jamais ecrase)."""
        ani_path = filedialog.askopenfilename(
            title="Select the legacy .ANI file",
            filetypes=[("Legacy .ANI", "*.ANI *.ani"), ("All files", "*.*")])
        if not ani_path:
            return
        base, _ext = os.path.splitext(ani_path)
        default_prmag_path = base + ".prmag"
        prmag_path = filedialog.asksaveasfilename(
            title="Save as .prmag", initialfile=os.path.basename(default_prmag_path),
            initialdir=os.path.dirname(default_prmag_path),
            defaultextension=".prmag", filetypes=[("STARpaleomag_Py .prmag", "*.prmag"), ("All files", "*.*")])
        if not prmag_path:
            return
        try:
            n_prmag = create_prmag_from_legacy_ani(ani_path, prmag_path)
        except OSError as e:
            self._showerror("Error", f"Could not read {ani_path}:\n{e}")
            return
        msg = (
            f"{n_prmag} specimen(s) -> {prmag_path}\n"
            "Only specimen id and core azimuth/dip/bedding come from the "
            ".ANI - site, date, geology, volume/mass are left as 'n.d'/"
            "defaults (fill in later with STARpaleomag_Py's Complete "
            "sample information...).\n"
        )

        pmagani_path = ani_path_for(prmag_path)
        if os.path.exists(pmagani_path):
            msg += f"companion .pmagani already exists, left untouched: {pmagani_path}\n"
        else:
            import_legacy_ani(ani_path, pmagani_path)
            msg += f"Also converted the tensors themselves -> {pmagani_path}\n"

        self.prmag_specimens = read_prmag_specimens(prmag_path)
        self.prmag_path = prmag_path
        self.donnees = read_ani_file(pmagani_path)
        self.ani_path = pmagani_path
        self.selection = []
        self._showinfo("prmag created from legacy .ANI", msg)
        self._afficher(msg)

    def ouvrir_creer_prmag_from_asc_dialog(self):
        """Cree un .prmag ET son .pmagani compagnon directement a partir
        d'un rapport .asc AGICO Kappabridge - demande explicite
        utilisateur ("ajouter en dessous de create prmag & pmagani from
        legacy .ANI... : create prmag & pmagani from .asc AGICO file"),
        equivalent de ouvrir_creer_prmag_from_ani_dialog mais lisant le
        fichier de mesure natif du kappabridge plutot qu'un .ANI deja
        converti - voir ams_prmag.create_prmag_from_asc pour le detail
        (colonnes cin/caz/dip/str_, meme deduplication par specimen
        qu'un .ANI : un specimen ATRM/AARM repete une ligne par variante
        jackknife, toutes partageant la meme orientation).

        Contrairement a la version .ANI (qui laisse un .pmagani DEJA
        existant intact plutot que d'ecraser), le compagnon est ici
        TOUJOURS construit via archive_asc_file : celle-ci est deja
        incrementale par construction (n'ajoute que les specimens
        vraiment nouveaux, ne touche jamais un .pmagani deja enrichi
        cote AMS_Py - moyennes de site, colonne export) - reproduire ici
        le garde-fou de la version .ANI serait redondant et empecherait
        un reimport legitime (nouvelles semaines de mesure sur le meme
        .asc, voir ouvrir_archiver_asc_dialog)."""
        asc_path = filedialog.askopenfilename(
            title="Select the AGICO .asc file",
            filetypes=[("AGICO .asc", "*.asc *.ASC"), ("All files", "*.*")])
        if not asc_path:
            return
        base, _ext = os.path.splitext(asc_path)
        default_prmag_path = base + ".prmag"
        prmag_path = filedialog.asksaveasfilename(
            title="Save as .prmag", initialfile=os.path.basename(default_prmag_path),
            initialdir=os.path.dirname(default_prmag_path),
            defaultextension=".prmag", filetypes=[("STARpaleomag_Py .prmag", "*.prmag"), ("All files", "*.*")])
        if not prmag_path:
            return
        try:
            n_prmag, warnings = create_prmag_from_asc(asc_path, prmag_path)
        except OSError as e:
            self._showerror("Error", f"Could not read {asc_path}:\n{e}")
            return
        msg = (
            f"{n_prmag} specimen(s) -> {prmag_path}\n"
            "Only specimen id and core azimuth/dip/bedding come from the "
            ".asc - site, date, geology, volume/mass are left as 'n.d'/"
            "defaults (fill in later with STARpaleomag_Py's Complete "
            "sample information...).\n"
        )
        if warnings:
            msg += f"{len(warnings)} block(s) skipped:\n" + "\n".join(f"  {w}" for w in warnings[:20])
            if len(warnings) > 20:
                msg += f"\n  ... and {len(warnings) - 20} more"
            msg += "\n"

        pmagani_path = ani_path_for(prmag_path)
        # `warnings` ci-dessus deja rapportes (meme fichier, meme parse
        # que create_prmag_from_asc) - archive_asc_file reparse le meme
        # .asc en interne (necessaire : elle a besoin des tenseurs, pas
        # seulement de l'orientation) mais ses PROPRES warnings ne sont
        # pas re-affiches, ce serait un doublon.
        pmagani_path, new_measurements, already_present, _asc_warnings = archive_asc_file(asc_path, pmagani_path)
        msg += (
            f"Tensors archived -> {pmagani_path}\n"
            f"{len(new_measurements)} new specimen(s), {len(already_present)} already present (skipped).\n"
        )

        self.prmag_specimens = read_prmag_specimens(prmag_path)
        self.prmag_path = prmag_path
        self.donnees = read_ani_file(pmagani_path)
        self.ani_path = pmagani_path
        self.selection = []
        self._showinfo("prmag created from .asc", msg)
        self._afficher(msg)

    def ouvrir_archiver_asc_dialog(self):
        """Archive un fichier .asc AGICO (rapport texte du kappabridge,
        logiciel SUSAR/Anisoft) dans un .pmagani - port de `readasc` (voir
        ams_asc.archive_asc_file pour le detail/l'etat de verification
        contre de vrais fichiers .asc) - demande explicite utilisateur
        ("can you also change the import asc files to new format", puis
        "l'acquisition de donnees et son archivage lors de la mesure dans
        le .asc est progressive et peut se faire sur plusieurs
        semaines... changer en archiver asc dans pmagani") : contrairement
        a l'ancienne "convert ASC to .pmagani" (qui REECRIVAIT tout le
        fichier a chaque passage), ce dialogue peut cibler un .pmagani
        DEJA existant (deja mesure les semaines precedentes, deja
        eventuellement complete cote AMS_Py - moyennes de site, colonne
        export) : seuls les specimens du .asc PAS ENCORE dans ce fichier
        sont ajoutes, le reste du fichier reste intact. cin/caz/dip/str_
        et les statistiques de Hext (err.%/F/F12/F23, quand presentes dans
        le rapport) sont ecrites directement (elles viennent reellement du
        .asc, contrairement au cas .prmag)."""
        path = filedialog.askopenfilename(
            title="Import .asc", filetypes=[("AGICO .asc", "*.asc *.ASC"), ("All files", "*.*")])
        if not path:
            return
        # Si un .pmagani est DEJA ouvert, archiver directement dedans -
        # demande explicite utilisateur ("si on a deja ouvert le
        # pmagani. Ensuite si on demande l'archivage, il suffit de
        # verifier que le meme specimen est bien present dans le prmag
        # associe") : pas besoin de redemander une cible, c'est deja
        # celle sur laquelle l'utilisateur travaille - la verification
        # utile se fait plus loin (specimens nouvellement archives
        # absents du .prmag compagnon), pas au choix du fichier.
        if self.ani_path and os.path.exists(self.ani_path):
            new_path = self.ani_path
        else:
            default_new_path = ani_path_for(path)
            # "Open" (pas "Save As") pour cibler un .pmagani DEJA
            # existant - demande explicite utilisateur ("peut on
            # archiver les nouvelles donnees d'un point .asc dans un
            # fichier existant sans avoir a le remplacer") :
            # asksaveasfilename affiche l'avertissement systeme macOS
            # "ce fichier existe deja, le remplacer ?" des qu'on choisit
            # un fichier existant - trompeur ici puisqu'on ne remplace
            # jamais rien (archive_asc_file n'ajoute que les specimens
            # manquants). Annuler cette 1ere boite -> propose d'en
            # creer un nouveau (asksaveasfilename, ou l'avertissement
            # est legitime : un NOUVEAU fichier qui porterait un nom
            # deja pris serait bien remplace).
            new_path = filedialog.askopenfilename(
                title="Archive into an existing .pmagani (Cancel to create a new file instead)",
                initialdir=os.path.dirname(default_new_path),
                filetypes=[("pmagani", "*.pmagani"), ("All files", "*.*")])
            if not new_path:
                new_path = filedialog.asksaveasfilename(
                    title="Create a new .pmagani",
                    initialfile=os.path.basename(default_new_path),
                    initialdir=os.path.dirname(default_new_path),
                    defaultextension=".pmagani", filetypes=[("pmagani", "*.pmagani"), ("All files", "*.*")])
                if not new_path:
                    return
        try:
            new_path, new_measurements, already_present, warnings = archive_asc_file(path, new_path)
        except OSError as e:
            self._showerror("Archive failed", f"Could not import {path}:\n{e}")
            return
        msg = (f"Archived {path} -> {new_path}\n"
               f"{len(new_measurements)} new specimen(s) archived, "
               f"{len(already_present)} already present (skipped).\n")
        if warnings:
            msg += f"{len(warnings)} block(s) skipped:\n" + "\n".join(f"  {w}" for w in warnings[:20])
            if len(warnings) > 20:
                msg += f"\n  ... and {len(warnings) - 20} more"
            msg += "\n"
        # self.donnees doit refleter l'etat COMPLET du .pmagani
        # maintenant que l'archivage est incremental (pas seulement les
        # measurements du .asc de cette session) - recharge tout depuis
        # le fichier, MAIS remplace chaque specimen present dans CE .asc
        # (nouveau ou deja_present) par sa version FRAICHEMENT PARSEE :
        # .pmagani n'a pas de colonnes cin/caz/dip/str_ (orientation),
        # donc une relecture pure perdrait celle des specimens que ce
        # .asc vient justement de fournir - le fichier lui-meme reste
        # inchange, seul self.donnees (etat en memoire de cette session)
        # beneficie de cette orientation le temps de cette session.
        asc_by_key = {
            (m.id.strip().upper(), m.code2.strip().upper()): m
            for m in new_measurements + already_present
        }
        self.donnees = [
            asc_by_key.get((m.id.strip().upper(), m.code2.strip().upper()), m)
            for m in read_ani_file(new_path)
        ]
        self.ani_path = new_path
        self.selection = []
        # orientation (cin/caz/dip/str) vient deja du .asc lui-meme ici -
        # PAS ecrasee par un .prmag compagnon (contrairement a ouvrir_ani_
        # dialog/ouvrir_import_legacy_ani_dialog) ; seul site/sample est
        # recupere du .prmag, pour "Export to Magic".
        prmag_candidate = self._prmag_path_for(new_path)
        if os.path.exists(prmag_candidate):
            self.prmag_specimens = read_prmag_specimens(prmag_candidate)
            self.prmag_path = prmag_candidate
            msg += f"companion prmag file found: {prmag_candidate}\n"
        # Verifie que les specimens NOUVELLEMENT archives figurent bien
        # dans le .prmag compagnon - demande explicite utilisateur (voir
        # plus haut) : un .pmagani est cense decrire des specimens du
        # .prmag associe (meme principe que calcul.find_orphan_ani_
        # specimens cote STARpaleomag_Py) ; un id absent du .prmag signale
        # un vrai probleme (faute de frappe, .asc d'un autre jeu de
        # donnees) plutot qu'un simple "pas encore mesure" - seuls les
        # specimens VRAIMENT ajoutes cette fois sont verifies (les
        # already_present l'ont deja ete lors d'un archivage precedent).
        if self.prmag_specimens and new_measurements:
            orphan_ids = sorted({
                m.id.strip() for m in new_measurements
                if self.prmag_specimens.get(m.id) is None
            })
            if orphan_ids:
                shown = ", ".join(orphan_ids[:20])
                more = f", ... ({len(orphan_ids) - 20} more)" if len(orphan_ids) > 20 else ""
                msg += (f"WARNING: {len(orphan_ids)} newly archived specimen(s) have NO "
                        f"counterpart in the companion .prmag - check for a mismatched "
                        f"file or a typo: {shown}{more}\n")
        self._afficher(msg)

    def lister_fichier_ani(self):
        """Equivalent de `listeani` (lect_asc.f:677-699)."""
        if not self.donnees:
            self._showwarning("No data", "Open a .pmagani/.ANI file first.")
            return
        import io
        buf = io.StringIO()
        list_measurements(self.donnees, out=buf)
        self._afficher(buf.getvalue())

    _ANI_MAGIC_INFO = {
        "A0": ("ATRM", "LP-AN-TRM"),
        "F0": ("AARM", "LP-AN-ARM"),
        "N0": ("AMS", "LP-X"),
    }

    # 1=echantillon(Sa)/2=in-situ(IS)/3=pendage corrige(TC) - MEME
    # etiquette que STARpaleomag_Py/calcul._ORIENT_MODE_TAG (pas un enum
    # CE/IS/CP invente separement - demande explicite utilisateur "CP
    # (use TC) and IS are in comment"), et meme convention que
    # self.orientation (radiobuttons "AMS data") / ams_selection.
    # _ORIENT_TO_FILE_CODE (le code REELLEMENT stocke dans la section
    # mean du .pmagani est le pourcentage MagIC 0-100, voir la-bas -
    # ce dict-ci n'est qu'un libelle d'AFFICHAGE console).
    _ORIENT_LABEL = {1: "Sa", 2: "IS", 3: "TC"}

    def exporter_magic_dialog(self):
        """Exporte les tenseurs 'A0'/'F0'/'N0' (le code MAJEUR de chaque
        experience - "with the A0 code being the major one" ; les variantes
        jackknife A+/A-/A1..B6 ne sont pas des enregistrements MagIC
        distincts) de self.donnees vers un fichier specimens.txt MagIC
        (aniso_s, colonne trace-normalisee a 1, format `s1:s2:s3:s4:s5:s6`
        - convention verifiee dans pmagpy.ipmag.calculate_aniso_parameters).
        Site/sample viennent de self.prmag_specimens (voir
        ouvrir_prmag_dialog) - demande explicite utilisateur ("it will be
        much easier to export Anisotropy data to Magic" une fois le .ANI
        associe a son .prmag).

        aniso_s_n_measurements/aniso_s_sigma/aniso_ftest/aniso_ftest12/
        aniso_ftest23/aniso_ftest_quality : simple passthrough des champs
        DEJA presents sur AMSMeasurement quand elle vient d'un .pmagani
        (n_positions/sigma/ftest/ftest12/ftest23/quality - statistiques
        de Hext calculees cote STARpaleomag_Py, voir ams_selection._read_pmagani_file)
        mais jusqu'ici jamais exportees - demande explicite utilisateur
        ("are there other parameters that were not previously calculated
        and useful for the magic export" -> "the easiest things to do").
        Vide (pas "n.d", convention MagIC pour une valeur manquante) pour
        une mesure qui ne les porte pas (ex. venant d'un ancien .ANI).

        aniso_tilt_correction = "-1" (coordonnees specimen) pour CHAQUE
        ligne : aniso_s est exporte tel que stocke, SANS reorientation
        (`self.orientation` n'est pas applique ici) - "-1" est d'ailleurs
        deja l'hypothese par defaut de pmagpy.ipmag lui-meme quand cette
        colonne est absente (`aniso_tilt_correction = -1  # assume
        specimen coordinates`), donc l'ecrire explicitement ne fait que
        rendre cette hypothese verifiable plutot que de laisser MagIC la
        deviner - colonne REQUISE des qu'un champ du groupe "Anisotropy"
        est renseigne (`requiredIfGroup`), absente de cet export jusqu'ici."""
        if not self.donnees:
            self._showwarning("No data", "Open a .pmagani/.ANI file first.")
            return
        rows = [m for m in self.donnees if m.code2 in self._ANI_MAGIC_INFO]
        if not rows:
            self._showwarning(
                "No data", "No A0/F0/N0 tensor found in the loaded .ANI data.")
            return
        path = filedialog.asksaveasfilename(
            title="Export to MagIC (specimens.txt)", defaultextension=".txt",
            initialfile="specimens.txt",
            filetypes=[("MagIC specimens", "*.txt"), ("All files", "*.*")])
        if not path:
            return

        header = [
            "specimen", "sample", "site", "aniso_type", "aniso_s", "aniso_tilt_correction",
            "aniso_s_n_measurements", "aniso_s_sigma",
            "aniso_ftest", "aniso_ftest12", "aniso_ftest23", "aniso_ftest_quality",
            "method_codes",
        ]
        lines = ["tab delimited\tspecimens", "\t".join(header)]
        missing_prmag = []
        exported = 0
        for m in rows:
            trace = m.k11 + m.k22 + m.k33
            if trace == 0:
                continue
            aniso_type, method_codes = self._ANI_MAGIC_INFO[m.code2]
            s1, s2, s3 = m.k11 / trace, m.k22 / trace, m.k33 / trace
            s4, s5, s6 = m.k12 / trace, m.k23 / trace, m.k13 / trace
            aniso_s = f"{s1:.6f}:{s2:.6f}:{s3:.6f}:{s4:.6f}:{s5:.6f}:{s6:.6f}"
            spec = self.prmag_specimens.get(m.id)
            if spec is None:
                missing_prmag.append(m.id)
            sample = spec.sample if spec else ""
            site = spec.site if spec else ""
            n_meas = str(m.n_positions) if m.n_positions is not None else ""
            sigma = f"{m.sigma:.6g}" if m.sigma is not None else ""
            ftest = f"{m.ftest:.6g}" if m.ftest is not None else ""
            ftest12 = f"{m.ftest12:.6g}" if m.ftest12 is not None else ""
            ftest23 = f"{m.ftest23:.6g}" if m.ftest23 is not None else ""
            quality = m.quality or ""
            lines.append("\t".join([
                m.id, sample, site, aniso_type, aniso_s, "-1",
                n_meas, sigma, ftest, ftest12, ftest23, quality,
                method_codes,
            ]))
            exported += 1

        with open(path, "w", encoding="utf-8", newline="") as f:
            f.write("\n".join(lines) + "\n")

        msg = f"MagIC specimens file written: {path}\n{exported} tensor(s) exported.\n"
        if missing_prmag:
            shown = ", ".join(missing_prmag[:10])
            more = "..." if len(missing_prmag) > 10 else ""
            msg += (f"No sample/site found for {len(missing_prmag)} specimen(s) - "
                    f"open the companion .prmag file first ({shown}{more}).\n")
        self._afficher(msg)

    def ouvrir_export_site_means_magic_dialog(self):
        """Exporte self.mean_results (tenseurs moyens DEJA calcules -
        "Tensorial mean..." - ou charges depuis un .pmagani via "Select
        results") vers un fichier sites.txt MagIC, au niveau SITE plutot
        que specimen (contrairement a "Export to Magic", qui exporte
        aniso_s specimen par specimen) - demande explicite utilisateur
        ("check in magic how to export a mean tensor at site level").
        Colonnes aniso_v1/v2/v3 (Matrix, un par axe propre), parametres
        de forme (aniso_p/l/f/t/pp/perc/total/ll/ff/vg/fl) et
        aniso_tilt_correction (voir ams_stats.MAGIC_TILT_CORRECTION_CODE,
        depuis l'orientation DEJA connue de chaque resultat - demande
        explicite utilisateur "are there other parameters that were not
        previously calculated and useful for the magic export" :
        colonne REQUISE par le data model des qu'un champ du groupe
        "Anisotropy" est renseigne, absente de cet export jusqu'ici) -
        voir ams_stats.magic_site_aniso_fields pour le detail exact des
        formules et l'ecart documente avec le formalisme Hext (1963) de
        pmagpy pour eta/zeta.

        Un resultat dont l'id ne se resout PAS a un site unique (voir
        ams_selection.mean_result_site - moyenne melangeant plusieurs
        specimens/sites) est ECARTE et journalise plutot qu'exporte sous
        un nom de site devine.

        `aniso_type`/`method_codes` (A0=ATRM/F0=AARM/N0=AMS, meme
        convention que `exporter_magic_dialog`) sont demandes UNE FOIS
        pour tout l'export, PAS par resultat : `TensorialMeanResult` ne
        porte pas le code2 du tenseur d'origine (`tsmean`/le format
        legacy "R ..." ne l'enregistrent pas non plus - voir
        ams_selection.read_ani_mean_results), donc rien ne permet de le
        retrouver automatiquement resultat par resultat."""
        if not self.mean_results:
            self._showwarning(
                "No mean result",
                "No mean tensor result in memory - run \"Tensorial mean...\" "
                "or \"Select results\" first.")
            return
        code_s = self._console_input(
            "Tensor type for all these means (A0=ATRM, F0=AARM, N0=AMS): ", "N0")
        if code_s is None:
            return
        code_s = code_s.strip().upper() or "N0"
        aniso_type, method_codes = self._ANI_MAGIC_INFO.get(code_s, ("AMS", "LP-X"))

        path = filedialog.asksaveasfilename(
            title="Export site means to MagIC (sites.txt)", defaultextension=".txt",
            initialfile="sites.txt",
            filetypes=[("MagIC sites", "*.txt"), ("All files", "*.*")])
        if not path:
            return

        header = [
            "site", "aniso_type", "aniso_tilt_correction", "aniso_v1", "aniso_v2", "aniso_v3",
            "aniso_p", "aniso_pp", "aniso_t", "aniso_l", "aniso_f",
            "aniso_perc", "aniso_total", "aniso_ll", "aniso_ff", "aniso_vg", "aniso_fl",
            "method_codes", "description",
        ]
        lines = ["tab delimited\tsites", "\t".join(header)]
        skipped = []
        exported = 0
        for orientation, res in self.mean_results:
            site = mean_result_site(res.id)
            if site is None:
                skipped.append(res.id or "(no id)")
                continue
            fields = magic_site_aniso_fields(res, aniso_type=aniso_type, orientation=orientation)
            row = [site] + [fields[key] for key in header[1:-2]] + [method_codes, fields["description"]]
            lines.append("\t".join(str(v) for v in row))
            exported += 1

        with open(path, "w", encoding="utf-8", newline="") as f:
            f.write("\n".join(lines) + "\n")

        msg = f"MagIC sites file written: {path}\n{exported} site mean(s) exported.\n"
        if skipped:
            shown = ", ".join(skipped[:10])
            more = "..." if len(skipped) > 10 else ""
            msg += (f"Skipped {len(skipped)} mean result(s) that mix more than one "
                    f"specimen/site prefix ({shown}{more}).\n")
        self._afficher(msg)

    def ouvrir_marquer_export_dialog(self):
        """Marque, dans le .pmagani couramment ouvert, la colonne "export"
        (voir ams_selection.mark_pmagani_export) : export=Y pour tout
        specimen present dans self.selection et toute moyenne presente
        dans self.mean_results, export=N pour TOUTE AUTRE ligne du meme
        fichier - demande explicite utilisateur ("is it possible to
        export from AMS_py only the data and mean tensors that we want
        to export and only these selected data will be taken into
        account in the main export from Starpaleomag"). self.selection/
        self.mean_results sont construits normalement (Select measurements/
        Tensorial mean/Select results) AVANT d'appeler ce menu - c'est cet
        ensemble courant qui devient le nouvel ensemble exporte, pas un
        ajout au precedent ("replace" et non "append", pour que ce menu
        reste rejouable a volonte sans devoir d'abord "deselectionner").
        STARpaleomag_Py lit ensuite cette meme colonne pour ignorer les
        lignes export=N lors de son propre "Export to Magic"."""
        if not self.ani_path or not os.path.isfile(self.ani_path):
            self._showwarning("No file", "Open a .pmagani file first.")
            return
        if not self.selection and not self.mean_results:
            self._showwarning(
                "Nothing selected",
                "Select measurements and/or compute mean results first - "
                "the current selection becomes the export set.")
            return

        selected_ids = {m.id.strip().upper() for m in self.selection if m.id}

        selected_mean_keys = set()
        skipped = []
        if self.mean_results:
            code_s = self._console_input(
                "Tensor type for all these means (A0=ATRM, F0=AARM, N0=AMS): ", "N0")
            if code_s is None:
                return
            default_code2 = code_s.strip().upper() or "N0"
            for orientation, res in self.mean_results:
                site = mean_result_site(res.id)
                if site is None:
                    skipped.append(res.id or "(no id)")
                    continue
                code2 = (res.source_code2 or default_code2).strip().upper()
                tilt_code = _ORIENT_TO_FILE_CODE.get(orientation, "")
                selected_mean_keys.add((site.upper(), code2, tilt_code))

        n_spec_in, n_spec_out, n_mean_in, n_mean_out = mark_pmagani_export(
            self.ani_path, selected_ids, selected_mean_keys)

        msg = (f"Export selection marked in {self.ani_path}\n"
               f"Specimens: {n_spec_in} kept (export=Y), {n_spec_out} excluded (export=N)\n"
               f"Means: {n_mean_in} kept (export=Y), {n_mean_out} excluded (export=N)\n")
        if skipped:
            shown = ", ".join(skipped[:10])
            more = "..." if len(skipped) > 10 else ""
            msg += (f"Skipped {len(skipped)} mean result(s) that mix more than one "
                    f"specimen/site prefix ({shown}{more}).\n")
        self._afficher(msg)

    # ------------------------------------------------------------------
    # AMS data
    # ------------------------------------------------------------------

    def selectionner_mesures(self):
        """Equivalent de `selmes` (lect_asc.f:436-561). ACCUMULE sur la
        selection existante - `nbmes` n'est JAMAIS remis a zero dans
        `selmes` lui-meme cote Fortran (seul `initmes`/"init list to
        zero" le fait, confirme en relisant le source : `nbmes=nbmes+1`
        pour chaque mesure trouvee, jamais `nbmes=0`) - meme bug/
        correction que cote STARpaleomag_Py ("the selection of results should
        not initialize the previous list"). Permet en particulier de
        constituer une selection couvrant les 15 lignes d'un specimen en
        plusieurs appels (ex. une plage d'etapes a la fois) sans perdre
        ce qui a deja ete selectionne, ou de cumuler plusieurs
        specimens avant un calcul groupe (Tensorial mean...) - demande
        explicite utilisateur ("select samples which permit especially
        to select the 15 lines per samples")."""
        if not self.donnees:
            self._showwarning("No data", "Open a .pmagani/.ANI file first (AMS Files > Open File .pmagani).")
            return
        pattern = self._console_input("Sample number: ", "*")
        if pattern is None:
            return
        smin = self._console_input("Step Min: ", "0")
        if smin is None:
            return
        smax = self._console_input("Step Max: ", "9000")
        if smax is None:
            return
        code = self._console_input("Code NRM Th or AF (example N0): ", "*")
        if code is None:
            return
        code = code.strip()
        demag1, demag2 = ("*", "*") if not code or code == "*" else (
            (code[0], "*") if len(code) == 1 else (code[0], code[1]))
        try:
            step_min = int(smin or 0)
            step_max = int(smax or 9000)
        except ValueError:
            self._showerror("Error", "Step Min/Max must be integers.")
            return
        new_matches = select_measurements(
            self.donnees, pattern or "*", step_min=step_min, step_max=step_max,
            demag1=demag1, demag2=demag2)
        self.selection = self.selection + new_matches
        self._afficher(f"+{len(new_matches)} data - total {len(self.selection)} data\n")
        self.lister_mesures()

    def lister_mesures(self):
        """Equivalent de `lismes`/`liste` (lect_asc.f:220-226,
        anisotropie.f:1250-1403) : PAS le tenseur brut (ca c'est `Liste File
        .ANI`/`listeani`) mais les axes propres k1/k2/k3 (dec/inc) + L/F/P/P'/T,
        calcules dans l'orientation couramment choisie."""
        if not self.selection:
            self._showwarning("No selection", "Select measurements first.")
            return
        self._afficher(format_measurement_list(self.selection, self.orientation.get()))

    def reinitialiser_liste(self):
        """Equivalent de `initmes` (AMS_OSX_x.f95:205-218)."""
        self.selection = []
        self._afficher("----\ndata list empty\n----\n")

    def supprimer_lignes(self):
        """Equivalent de `delmes` (lect_asc.f:563-590), UN SEUL numero a la
        fois a l'origine - etendu a une liste d'indices ("1,3,5-7", meme
        syntaxe que STARpaleomag_Py.ouvrir_delete_results_dialog) pour
        retirer plusieurs mesures en un seul passage - demande explicite
        utilisateur ("to delete measurements, perhaps best to uses the
        specimen numbers as in STARpaleomag_Py")."""
        if not self.selection:
            self._showwarning("No selection", "Select measurements first.")
            return
        import io
        buf = io.StringIO()
        for i, m in enumerate(self.selection, start=1):
            buf.write(f"{i:4d} : {m.id:<12}")
            if i % 3 == 0:
                buf.write("\n")
        self._afficher(buf.getvalue())
        indices_s = self._console_input(
            "Indices to delete (e.g. 1,3,5-7 - see the numbers above; "
            "Escape to cancel): ", "")
        if indices_s is None:
            return
        indices_s = indices_s.strip()
        if not indices_s:
            return

        to_delete = set()
        try:
            for token in indices_s.split(","):
                token = token.strip()
                if not token:
                    continue
                if "-" in token:
                    a, b = token.split("-", 1)
                    to_delete.update(range(int(a), int(b) + 1))
                else:
                    to_delete.add(int(token))
        except ValueError:
            self._showerror("Error", "Invalid index list (use e.g. 1,3,5-7).")
            return

        n = len(self.selection)
        invalid = sorted(i for i in to_delete if i < 1 or i > n)
        if invalid:
            self._showerror("Error", f"Index out of range: {invalid} (1-{n}).")
            return

        self.selection = [m for i, m in enumerate(self.selection, start=1) if i not in to_delete]
        self._afficher(f"{len(to_delete)} measurement(s) removed - {len(self.selection)} remaining.\n")
        self.lister_mesures()

    def changer_orientation(self):
        """Equivalent de `selech`/`selis`/`selcp` (AMS_OSX_x.f95:309-376)."""
        if self.selection:
            self.lister_mesures()
        self._refresh_current_graphic()

    # ------------------------------------------------------------------
    # Graphics
    # ------------------------------------------------------------------

    def _refresh_current_graphic(self):
        if self._current_graphic is None:
            return
        self.fig.clear()
        orientation = self.orientation.get()
        if self._current_graphic == "stereo":
            measurements = self.selection if self._aniso_source in (0, 2) else []
            means = ([r for _o, r in self.mean_results] if self._aniso_source in (1, 2) else [])
            if not measurements and not means:
                return
            self.fig.set_size_inches(5.5, 5.5, forward=True)
            build_stereo_figure(
                measurements, orientation=orientation, invert_negative=self.invert_negative,
                ams_iproj=self.ams_iproj, mean_results=means, fig=self.fig)
        elif self._current_graphic == "aniso_params":
            entries = []
            if self._aniso_source in (0, 2):
                entries += shape_entries_from_measurements(self.selection, orientation)
            if self._aniso_source in (1, 2):
                entries += shape_entries_from_mean_results([r for _o, r in self.mean_results])
            if not entries:
                return
            self.fig.set_size_inches(5.5, 9.5, forward=True)
            build_anisotropy_parameters_figure(
                entries, fig=self.fig,
                max_scale_top=self._aniso_scale_top, max_scale_bottom=self._aniso_scale_bottom)
        elif self._current_graphic == "susc_aniso":
            if not self.selection:
                return
            self.fig.set_size_inches(5.5, 9.5, forward=True)
            build_susceptibility_anisotropy_figure(self.selection, orientation, fig=self.fig)
        elif self._current_graphic == "im_re_susc":
            if not self.selection:
                return
            self.fig.set_size_inches(7.0, 5.5, forward=True)
            build_im_re_susceptibility_figure(self.selection, fig=self.fig)
        elif self._current_graphic == "bootstrap":
            if self._last_bootstrap is None:
                return
            self.fig.set_size_inches(5.5, 5.5, forward=True)
            build_bootstrap_stereo_figure(
                self._last_bootstrap, orientation=orientation, invert_negative=self.invert_negative,
                ams_iproj=self.ams_iproj, show_cloud=self._bootstrap_show_cloud,
                show_ellipse=self._bootstrap_show_ellipse, fig=self.fig)
        else:
            return
        self._redraw_canvas()

    def _on_plot_pick(self, event):
        """Gestionnaire unique de clic-sur-artiste (pick_event) pour le
        graphique "Anisotropy parameters" - voir ams_xy._register_point_
        pick/_register_xaxis_pick pour ce qui est enregistre et pourquoi
        (meme mecanisme que STARpaleomag_Py app._on_plot_pick/xygraph.
        _pick_points). `event.ind` : toujours un seul point par artiste
        pour "aniso_specimen" (un scatter par point), d'ou l'usage de
        `event.ind[0]` - absent pour "aniso_xaxis" (l'artiste est l'axe
        lui-meme, pas un scatter). Le panneau (haut/bas) est identifie
        par la POSITION de l'Axes dans la figure (`self.fig.axes[0]` =
        add_subplot(211) = haut), pas par son contenu (Flinn ou P'/T-Im
        selon le cas - voir ams_xy.build_anisotropy_parameters_figure) :
        chaque panneau garde ainsi sa PROPRE echelle, correction demande
        explicite utilisateur ("the scale should not apply to both plots
        (Re et Im) as they have very different ranges")."""
        artist = event.artist
        kind = getattr(artist, "_ams_pick_kind", None)
        if kind == "aniso_specimen":
            data = getattr(artist, "_ams_pick_data", None)
            if data is not None and len(event.ind) > 0:
                self._afficher(f"[Anisotropy parameters] {data}\n")
        elif kind == "aniso_xaxis":
            is_top = artist.axes is self.fig.axes[0] if self.fig.axes else True
            self._prompt_single_aniso_scale("top" if is_top else "bottom")

    def _prompt_single_aniso_scale(
        self, panel: str, refresh: bool = True, reset_default: bool = False,
    ) -> bool:
        """Invite pour changer la borne haute d'echelle (F/L ou P') d'UN
        SEUL panneau ("top"/"bottom") de "Anisotropy parameters" - appelee
        soit deux fois de suite a l'ouverture du graphique (voir
        afficher_anisotropy_parameters, `refresh=False` : c'est l'appelant
        qui redessine une fois les deux valeurs saisies), soit une seule
        fois depuis un clic sur l'axe X d'un panneau (`refresh=True`,
        demande explicite utilisateur "is it possible to click on the X
        axis to change the max value"). Retourne False si l'utilisateur
        annule (Escape) ou entre une valeur invalide - l'appelant ne doit
        alors pas continuer.

        `reset_default` (utilise par afficher_anisotropy_parameters) :
        n'affiche PAS la valeur precedente comme suggestion, meme si elle
        est encore memorisee - demande explicite utilisateur ("une fois
        que le choix de l'echelle a ete fait, il n'est pas remis a zero,
        contrairement a la question donnees(0) resultats(1) d+r(2)") :
        cette derniere propose TOUJOURS "0" par defaut a l'ouverture,
        jamais le dernier choix - l'echelle doit se comporter pareil a
        l'ouverture du graphique. Un clic sur l'axe X (refresh=True), en
        revanche, ajuste une valeur EXISTANTE : y montrer le dernier
        choix comme suggestion reste logique et n'est pas concerne."""
        attr = "_aniso_scale_top" if panel == "top" else "_aniso_scale_bottom"
        current = None if reset_default else getattr(self, attr)
        scale_default = "" if current is None else f"{current:g}"
        scale_s = self._console_input(
            f"Max scale for {panel} panel (P>{P_OUTLIER_THRESHOLD:g} always excluded as outliers, "
            "empty = auto) : ", scale_default)
        if scale_s is None:
            return False
        scale_s = scale_s.strip()
        if not scale_s:
            setattr(self, attr, None)
        else:
            try:
                setattr(self, attr, float(scale_s))
            except ValueError:
                self._showerror("Error", "Scale must be a number.")
                return False
        if refresh:
            self._refresh_current_graphic()
        return True

    def _redraw_canvas(self):
        """Port de STARpaleomag_Py/app._redraw_canvas - MEME probleme rencontre
        ici ("in AMS_Py, we have the same pb encountered with starmac
        about the clear screen and the initial position of the graphs,
        corrected by moving the window" - demande explicite utilisateur) :
        ni `update_idletasks()` ni `sashpos()`/`geometry()` seuls ne
        forcent un vrai cycle <Configure> sur le widget canvas Tk - seul
        un VRAI redimensionnement le fait, d'ou le "secouement" de la
        fenetre (+1px puis retour) qui remplace le geste manuel de
        deplacer/redimensionner la fenetre que l'utilisateur faisait pour
        forcer l'affichage. Agrandit le volet gauche du PanedWindow/la
        fenetre pour accueillir la Figure demandee (plafonne a l'espace
        ecran disponible), secoue la fenetre, puis aligne la Figure sur la
        taille REELLE obtenue par le widget canvas (source de verite,
        matplotlib resynchronise de toute facon dessus au moindre
        <Configure> ulterieur - inutile de lutter contre ca)."""
        width_in, height_in = self.fig.get_size_inches()
        dpi = self.fig.dpi

        max_width_px = max(400, self.root.winfo_screenwidth() - 100)
        max_height_px = max(300, self.root.winfo_screenheight() - 100)
        avail_fig_w_in = (max_width_px - 20) / dpi
        avail_fig_h_in = (max_height_px - 90) / dpi
        if width_in > avail_fig_w_in or height_in > avail_fig_h_in:
            scale = min(avail_fig_w_in / width_in, avail_fig_h_in / height_in, 1.0)
            width_in, height_in = width_in * scale, height_in * scale

        width_px = int(width_in * dpi) + 20
        height_px = int(height_in * dpi) + 90  # marge titre/menus/barre d'etat

        try:
            current_w = self.paned_window.sashpos(0)
            if width_px > current_w:
                self.paned_window.sashpos(0, width_px)
        except tk.TclError:
            pass

        if height_px > self.root.winfo_height():
            self.root.geometry(f"{self.root.winfo_width()}x{height_px}")

        self.root.update_idletasks()

        w = self.root.winfo_width()
        h = self.root.winfo_height()
        self.root.geometry(f"{w + 1}x{h + 1}")
        self.root.update_idletasks()
        self.root.geometry(f"{w}x{h}")
        self.root.update_idletasks()

        canvas_w = self.canvas_fig.get_tk_widget().winfo_width()
        canvas_h = self.canvas_fig.get_tk_widget().winfo_height()
        if canvas_w > 1 and canvas_h > 1:
            self.fig.set_size_inches(canvas_w / dpi, canvas_h / dpi, forward=False)

        self.canvas_fig.draw()

    def exporter_svg(self):
        """Equivalent de STARpaleomag_Py/app.exporter_svg (menu Graphics >
        Export SVG...) - demande explicite utilisateur ("Could you add the
        svg export in AMS_Py"), remplace le stub "save SVG graphics"
        (_not_implemented) deja present dans le menu. Meme approche que
        cote STARpaleomag_Py : pas besoin d'un code d'export separe,
        matplotlib ecrit un SVG directement depuis la Figure affichee
        (self.fig) - svgwriter.py (module partage entre les deux
        applications) sert a un pipeline de rendu different, independant
        de ce menu."""
        if not self.fig.axes:
            self._showwarning("No graphic", "Display a graphic first (e.g. Stereographique).")
            return

        default_name = f"{self._current_graphic}.svg" if self._current_graphic else "graphic.svg"
        path = filedialog.asksaveasfilename(
            title="Export as SVG",
            defaultextension=".svg",
            initialfile=default_name,
            filetypes=[("SVG", "*.svg"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            self.fig.savefig(path, format="svg")
        except Exception as e:
            self._showerror("Error", f"SVG export failed:\n{e}")
            return
        self._showinfo("Export successful", f"Graphic exported:\n{path}")

    def ouvrir_user_guide(self):
        """Ouvre le guide utilisateur (AMS_Py Guide) dans le navigateur
        systeme - equivalent de STARpaleomag_Py/app.ouvrir_user_guide (meme
        principe : un fichier HTML STATIQUE livre EN LOCAL avec l'appli,
        voir _resource_path, help/AMS_Py_Guide.html) - demande explicite
        utilisateur ("is it possible to write a guide for AMS_Py").

        Saute directement au bloc du DERNIER menu de contenu utilise
        (self._help_anchor, voir _menu_cmd) - ouvre la page d'accueil du
        guide (aucune ancre) tant qu'aucun item de menu de contenu n'a
        encore ete invoque dans cette session."""
        guide_path = _resource_path("help", "AMS_Py_Guide.html")
        if not os.path.exists(guide_path):
            self._showerror("Error", f"User guide not found:\n{guide_path}")
            return
        url = f"file://{guide_path}"
        if self._help_anchor:
            url += f"#{self._help_anchor}"
        webbrowser.open(url)

    def afficher_stereo(self):
        """Equivalent de `sterams`/`stereo`+`tratm` (menu Graphics >
        Stereographique) : axes k1/k2/k3 des mesures brutes ET/OU des
        resultats de moyenne tensorielle deja calcules (equivalent partiel
        de `tratm`, sans les ellipses de confiance - symboles 2x plus gros
        pour les distinguer). Meme prompt "donnees(0) resultats(1) d+r(2)"
        que `afficher_anisotropy_parameters` (etat partage, `self.
        _aniso_source`)."""
        if not self.selection and not self.mean_results:
            self._showwarning(
                "No data", "Select measurements and/or compute a tensorial mean first.")
            return
        src_s = self._console_input("data(0) results(1) d+r(2) : ", str(self._aniso_source))
        if src_s is None:
            return
        try:
            self._aniso_source = int(src_s)
        except ValueError:
            self._aniso_source = 0
        if self._aniso_source not in (0, 1, 2):
            self._aniso_source = 0
        self._current_graphic = "stereo"
        self._refresh_current_graphic()

    def afficher_anisotropy_parameters(self):
        """Equivalent (matplotlib, pas un port pixel-pres) de `flint`/`flinn`
        + `tpprim` (menu Graphics > Anisotropy parameters - entree UNIQUE,
        remplace les deux anciennes entrees Flinn/T-Pprim) : Flinn et P'/T
        affiches ensemble sur la meme page (meme motif que xygraph.py cote
        STARpaleomag_Py). Source des points : donnees brutes (0), resultats de
        moyenne tensorielle deja calcules (1), ou les deux (2) - meme
        prompt que le Fortran ("donnees(0) resultats(1) d+r(2)")."""
        if not self.selection and not self.mean_results:
            self._showwarning(
                "No data", "Select measurements and/or compute a tensorial mean first.")
            return
        src_s = self._console_input("data(0) results(1) d+r(2) : ", "0")
        if src_s is None:
            return
        try:
            self._aniso_source = int(src_s)
        except ValueError:
            self._aniso_source = 0
        if self._aniso_source not in (0, 1, 2):
            self._aniso_source = 0
        if not self._prompt_single_aniso_scale("top", refresh=False, reset_default=True):
            return
        if not self._prompt_single_aniso_scale("bottom", refresh=False, reset_default=True):
            return
        self._current_graphic = "aniso_params"
        self._refresh_current_graphic()

    def afficher_susceptibility_anisotropy(self):
        """Nouveau graphique (pas de source Fortran - demande explicite
        utilisateur "can we built a plot with two graphs of susceptibility
        on X and anisotropy degree on Y, above the Im and below the RE")
        - voir ams_xy.build_susceptibility_anisotropy_figure. Un point
        par mesure BRUTE de self.selection (pas de choix "donnees/
        resultats" : un tenseur moyen n'a pas de susceptibilite bulk
        individuelle comparable). Panneau du bas (Real) : susceptibilite
        en echelle log normale, min a gauche et max a droite ("from min
        to max"). Panneau du haut (Imaginary) :
        |susceptibilite negative| en echelle log, axe INVERSE - grandeur
        croissante VERS LA GAUCHE (demande explicite utilisateur "for the
        negative K, plotting with the scale increasing to the left").
        Y = P' (degre d'anisotropie corrige de Jelinek, meme convention
        que le panneau P'/T de "Anisotropy parameters")."""
        if not self.selection:
            self._showwarning("No data", "Select measurements first.")
            return
        self._current_graphic = "susc_aniso"
        self._refresh_current_graphic()

    def afficher_im_re_susceptibility(self):
        """Nouveau graphique (pas de source Fortran - demande explicite
        utilisateur "a final plot with the Im susceptibility on X and the
        Re susceptibility on Y same conventions as above for Kim + and
        -") - voir ams_xy.build_im_re_susceptibility_figure. Un point par
        SPECIMEN de self.selection ayant a la fois une mesure Real et une
        mesure Imaginary (susceptibilite scalaire, independante de
        l'orientation - pas de prompt d'orientation ici). Memes
        conventions que "Susceptibility / anisotropy degree" pour le
        cote Im : K negatif a gauche (|susceptibilite| en log, axe
        inverse), K positif a droite (susceptibilite en log, axe normal),
        echelles X independantes, ligne verticale a leur jonction - Y
        (susceptibilite Real) egalement en log, partagee entre les deux
        moities."""
        if not self.selection:
            self._showwarning("No data", "Select measurements first.")
            return
        self._current_graphic = "im_re_susc"
        self._refresh_current_graphic()

    def ouvrir_stereo_sample_dialog(self):
        """Equivalent de `sterams2` (AMS_OSX_x.f95:239-296, menu Graphics >
        Stereographic sample) : parcourt self.selection SPECIMEN PAR
        SPECIMEN (regroupement par id consecutif, meme logique que
        format_measurement_list) et affiche pour chacun, dans l'ordre, sa
        liste (avec le commentaire A0), son stereogramme, puis son
        diagramme de forme (Flinn/P'-T) - avec les memes invites que le
        Fortran ("return for T-Pprim" puis "next sample: return or
        (q)uit"). Escape a n'importe quelle invite interrompt le
        parcours."""
        if not self.selection:
            self._showwarning("No data", "Select measurements first.")
            return
        groups = []
        for m in self.selection:
            if groups and groups[-1][0] == m.id:
                groups[-1][1].append(m)
            else:
                groups.append((m.id, [m]))

        self.text_area.insert(tk.END, "\n--- Stereographic sample (Escape to cancel) ---\n", "prompt")
        orientation = self.orientation.get()
        self._current_graphic = None  # ce parcours dessine lui-meme, sans passer par _refresh_current_graphic
        for idx, (specimen_id, group) in enumerate(groups):
            self._afficher(
                f"--- {specimen_id} ({idx + 1}/{len(groups)}) ---\n"
                + format_measurement_list(group, orientation)
            )
            self.fig.clear()
            self.fig.set_size_inches(5.5, 5.5, forward=True)
            build_stereo_figure(
                group, orientation=orientation, invert_negative=self.invert_negative,
                ams_iproj=self.ams_iproj, fig=self.fig)
            self._redraw_canvas()

            choice = self._console_input(f"[{specimen_id}] return for T-Pprim (Escape to stop): ", "")
            if choice is None:
                return

            entries = shape_entries_from_measurements(group, orientation)
            self.fig.clear()
            self.fig.set_size_inches(5.5, 9.5, forward=True)
            build_anisotropy_parameters_figure(
                entries, fig=self.fig,
                max_scale_top=self._aniso_scale_top, max_scale_bottom=self._aniso_scale_bottom)
            self._redraw_canvas()

            if idx == len(groups) - 1:
                break
            choice = self._console_input(f"[{specimen_id}] next sample: return or (q)uit: ", "")
            if choice is None or choice.strip().lower() == "q":
                break

    def ouvrir_stereo_site_dialog(self):
        """Meme parcours que `ouvrir_stereo_sample_dialog` ("Stereographic
        sample"), mais regroupe self.selection SITE PAR SITE au lieu de
        specimen par specimen - demande explicite utilisateur ("est ce
        possible d'ajouter l'option stereo sites pour une visualisation
        rapide des donnees par site, suivant l'idee de stereo sample").
        Site = les 6 premiers caracteres de l'id (annee(2)+site(4)) -
        meme convention que `ams_selection.mean_result_site`/
        STARpaleomag_Py field_notes.py (ex. "24WH0103A"/"24WH0104B" ->
        site "24WH01"), regroupement par PREFIXE CONSECUTIF (comme le
        regroupement par id de la version specimen). Un site typique
        contenant plusieurs specimens (chacun potentiellement Re+Im),
        le stereogramme et le Flinn/P'-T affiches par site font
        naturellement basculer ce dernier sur la disposition P'-Im/P'-Re
        separee des que le site melange les deux (voir ams_xy.
        build_anisotropy_parameters_figure) - utile precisement pour ce
        genre de coup d'oeil rapide multi-specimens."""
        if not self.selection:
            self._showwarning("No data", "Select measurements first.")
            return
        groups = []
        for m in self.selection:
            site = m.id[:6]
            if groups and groups[-1][0] == site:
                groups[-1][1].append(m)
            else:
                groups.append((site, [m]))

        self.text_area.insert(tk.END, "\n--- Stereographic site (Escape to cancel) ---\n", "prompt")
        orientation = self.orientation.get()
        self._current_graphic = None  # ce parcours dessine lui-meme, sans passer par _refresh_current_graphic
        for idx, (site, group) in enumerate(groups):
            self._afficher(
                f"--- {site} ({idx + 1}/{len(groups)}, {len(group)} measurement(s)) ---\n"
                + format_measurement_list(group, orientation)
            )
            self.fig.clear()
            self.fig.set_size_inches(5.5, 5.5, forward=True)
            build_stereo_figure(
                group, orientation=orientation, invert_negative=self.invert_negative,
                ams_iproj=self.ams_iproj, fig=self.fig)
            self._redraw_canvas()

            choice = self._console_input(f"[{site}] return for T-Pprim (Escape to stop): ", "")
            if choice is None:
                return

            entries = shape_entries_from_measurements(group, orientation)
            self.fig.clear()
            self.fig.set_size_inches(5.5, 9.5, forward=True)
            build_anisotropy_parameters_figure(
                entries, fig=self.fig,
                max_scale_top=self._aniso_scale_top, max_scale_bottom=self._aniso_scale_bottom)
            self._redraw_canvas()

            if idx == len(groups) - 1:
                break
            choice = self._console_input(f"[{site}] next site: return or (q)uit: ", "")
            if choice is None or choice.strip().lower() == "q":
                break

    def ouvrir_paramster_dialog(self):
        """Equivalent de `paramster` (anisotropie.f:2166-2208) : projection
        (stereo/lambert) et inversion des inclinaisons negatives pour le
        stereonet. Le rayon personnalise ("rayon de la projection") n'est
        PAS demande : le Fortran l'ignore de toute facon (bug reel confirme,
        `rdim` reste toujours 15.0 - voir docstring ams_stereo.py)."""
        proj_s = self._console_input(
            "projection stereo(1)  lambert(0) : ", str(self.ams_iproj))
        if proj_s is None:
            return
        try:
            self.ams_iproj = 1 if int(proj_s) == 1 else 0
        except ValueError:
            pass
        inv_s = self._console_input(
            "inversion des inclinaisons negatives ?(y/n): ", "y" if self.invert_negative else "n")
        if inv_s is None:
            return
        self.invert_negative = not inv_s.strip().lower().startswith("n")
        self._refresh_current_graphic()

    # ------------------------------------------------------------------
    # Calcul
    # ------------------------------------------------------------------

    def ouvrir_corpaldir_dialog(self):
        """Equivalent de `corpaldir` (Util_AMS.f95:2-115) : correction
        inverse d'une direction paleomagnetique par le PREMIER tenseur de la
        selection courante ("attention prend le premier tenseur de la
        liste")."""
        if not self.selection:
            self._showwarning("No selection", "Select measurements first (need at least 1 tensor).")
            return
        self._afficher(
            " correction of Pmag data with TRM tensor\n"
            " attention prend le premier tenseur de la liste\n"
            " attention de bien choisir en CE ou IS selon la direction a corriger\n"
        )
        vals = self._console_input("Intensity, Declinaison, Inclinaison : ", "")
        if vals is None:
            return
        try:
            intensity, dec, inc = (float(v) for v in vals.replace(",", " ").split())
        except ValueError:
            self._showerror("Error", "Enter 3 numbers: intensity dec inc")
            return
        tensor = self.selection[0]
        x1, y1, z1 = correct_direction_with_tensor(tensor, intensity, dec, inc)
        self._afficher(f"Inverse Correction Tensor: Int:{x1:5.1f} Dec:{y1:6.1f} Inc:{z1:6.1f}\n")

    def ouvrir_soustract_dialog(self):
        """Equivalent de `soustract` pour exactement 2 mesures (Util_AMS.f95:158-283)
        - la variante liste (`soustractliste`, plusieurs paires) n'est pas
        encore portee."""
        if len(self.selection) != 2:
            self._showwarning(
                "Invalid selection",
                "Select exactly 2 measurements (soustractliste, for several "
                "consecutive pairs, is not yet ported).")
            return
        pct_s = self._console_input("% subtraction between 1 and 100 : ", "100")
        if pct_s is None:
            return
        try:
            percent = float(pct_s)
        except ValueError:
            self._showerror("Error", "Must be a number.")
            return
        res = subtract_tensors(self.selection[0], self.selection[1], percent=percent)
        out = []
        for m in (res.diff_normalized, res.diff_inverse):
            out.append(
                f"D {m.id:<12} {m.cin:.4f} {m.caz:.4f} {m.dip:.4f} {m.str_:.4f} 1 {m.code2} "
                f"{m.k11:.5f} {m.k22:.5f} {m.k33:.5f} {m.k12:.5f} {m.k23:.5f} {m.k13:.5f} {m.s:.4f}"
            )
        self._afficher(
            "\n----------------\ndifference between two tensors\n----------------\n"
            + "\n".join(out) + "\n"
        )
        self.selection = [self.selection[0], self.selection[1], res.diff_weighted,
                           res.diff_normalized, res.diff_inverse]
        self.lister_mesures()

    def ouvrir_tsmean_dialog(self):
        """Equivalent de `tsmean` (anisotropie.f:1668-1801) : moyenne
        tensorielle avec statistiques de Jelinek (1978) sur la selection
        courante, dans l'orientation actuellement choisie (AMS data > Sample/
        In situ/Tilt corrected). Necessite au moins 3 mesures.

        VERIFIE : les parametres de forme (L,F,P,P') via `ams_stats.
        shape_params` sont confirmes exacts contre magicams.txt reel
        (specimen 98PL0401A). NON VERIFIE de bout en bout : le pipeline
        complet erbar/ellips (ellipses de confiance a 95%) - transcrit
        formule a formule depuis le source, mais sans jeu de donnees reel
        complet (mêmes specimens + orientation exacte) pour comparaison
        octet-pres. A traiter avec prudence pour les angles de confiance
        tant qu'un cas reel verifiable n'est pas disponible.

        Reel/Imaginaire separes (demande explicite utilisateur - "pour le
        calcul des tenseurs moyens, il vaut mieux separer Re et Im car
        les degres d'anisotropie sont tres differents, quand les deux
        sont selectionnes, calculer 2 tenseurs") : si `self.selection`
        melange au moins une mesure reelle et une imaginaire (voir
        `is_imaginary_component`), DEUX moyennes tensorielles sont
        calculees separement (Imaginary puis Real) plutot qu'une seule
        moyenne melangeant des degres d'anisotropie incomparables -
        chacune gardee en memoire independamment sur confirmation. Une
        selection homogene (tout Re, tout Im) garde le comportement a
        un seul tenseur, inchange.

        Correction kmax/kmin des tenseurs a susceptibilite negative
        AVANT moyenne (demande explicite utilisateur - "les tenseurs
        moyens de Im negatifs ont aussi le kmax et le kmin intervertis.
        Comment faire une moyenne pour un site avec une partie des
        echantillons avec des Kim negatifs et d'autres positifs") : voir
        `ams_stats.negative_k_corrected_tensor`, applique ici a chaque
        mesure ou `m.s<0` (tenseur entier deja retourne par Agico) avant
        de le donner a `tsmean` - necessaire des qu'un groupe mele des
        specimens negatifs et positifs (leur "grande valeur propre" ne
        represente pas le meme axe physique), pas seulement pour
        l'affichage individuel (`ams_stats.principal_axes`, qui ne fait
        que relabeliser APRES coup, insuffisant pour une moyenne).

        Outliers (P=k1/k3 > `ams_xy.P_OUTLIER_THRESHOLD`, ou indefini -
        voir `ams_stats.shape_params`) EXCLUS de la moyenne, PAS
        seulement des graphiques (demande explicite utilisateur - "in
        the mean tensor, the direction and ellipses are OK but the
        anisotropy degree just explose to irrealistic values" sur un
        site reel, 24WH36) : un seul specimen dont la partie imaginaire
        est dominee par le bruit (susceptibilite quasi-nulle, degre
        d'anisotropie individuel deja demesure - confirme jusque dans le
        propre calcul d'Agico sur le fichier .asc d'origine, PAS un bug
        d'analyse) suffit a faire tendre un axe de la moyenne vers zero,
        avec `normalize=True` (division par la trace propre de CE
        specimen, deja quasi-nulle) OU MEME sans normaliser (l'axe
        moyen lui-meme peut s'annuler par pure annulation directionnelle
        entre specimens faibles) - un ratio dont le denominateur est
        quasi-nul explose quelle que soit la normalisation choisie.
        Exclure ces specimens AVANT le calcul (verifie sur 24WH36 : P'
        chute de ~295-900 a 2.5, une valeur physiquement plausible) est
        la seule correction qui traite la cause reelle plutot que le
        symptome."""
        if len(self.selection) <= 2:
            self._showwarning(
                "Not enough data", "Tensorial mean needs at least 3 measurements (il>2, like the Fortran).")
            return
        orientation = self.orientation.get()
        if orientation == 1:
            confirm = self._console_input("Sample coordinate ; should we continue (Y/n)?: ", "Y")
            if confirm is None:
                return
            if confirm.strip().lower().startswith("n"):
                return
        norm_s = self._console_input("normalisation of the tensors (Y/n) ?: ", "Y")
        if norm_s is None:
            return
        normalize = not norm_s.strip().lower().startswith("n")

        im_group = [m for m in self.selection if is_imaginary_component(m)]
        re_group = [m for m in self.selection if not is_imaginary_component(m)]
        if im_group and re_group:
            groups = [("Imaginary", im_group), ("Real", re_group)]
            self._afficher(
                "Selection mixes Real and Imaginary tensors (very different anisotropy "
                "degrees) - computing 2 separate mean tensors.\n")
        else:
            groups = [(None, self.selection)]

        for label, group in groups:
            prefix = f"[{label}] " if label else ""
            filtered = []
            n_outliers = 0
            for m in group:
                axes = principal_axes(m, orientation)
                p = shape_params(axes[0][0], axes[1][0], axes[2][0])["P"]
                if math.isnan(p) or p > P_OUTLIER_THRESHOLD:
                    n_outliers += 1
                    continue
                filtered.append(m)
            if n_outliers:
                self._afficher(
                    f"{prefix}{n_outliers} outlier(s) excluded (P>{P_OUTLIER_THRESHOLD:g} or "
                    "undefined) before averaging.\n")
            if len(filtered) <= 2:
                self._afficher(
                    f"{prefix}skipped: needs at least 3 valid measurements (has {len(filtered)}).\n")
                continue
            tensors = []
            for m in filtered:
                a = apply_orientation(m, orientation)
                if m.s < 0.0:
                    a = negative_k_corrected_tensor(a)
                tensors.append((a[0, 0], a[1, 1], a[2, 2], a[0, 1], a[1, 2], a[0, 2]))
            res = tsmean(tensors, normalize=normalize)
            if res is None:
                if label is None:
                    self._showwarning("No result", "Mean tensor is isotropic or not enough valid tensors.")
                else:
                    self._afficher(f"{prefix}mean tensor is isotropic or not enough valid tensors - skipped.\n")
                continue
            # equivalent de `storeres` (anisotropie.f:2989-3049) :
            # id = 6 premiers caracteres du 1er specimen + 6 premiers du dernier
            res.id = filtered[0].id[:6] + filtered[-1].id[:6]
            if all(is_imaginary_component(m) for m in filtered):
                res.source_code2 = "IM"
            elif not any(is_imaginary_component(m) for m in filtered):
                res.source_code2 = "RE"
            if res.ellipsoid_type == 0:
                self._afficher(f"{prefix}Isotropic mean tensor - no axis statistics.\n")
                continue

            header = f"--- {label} ---\n" if label else ""
            self._afficher(header + format_tsmean_box(res, orientation))

            keep_s = self._console_input(
                f"Would you like to keep this {prefix}result in memory ?(Y/n): ", "Y")
            if keep_s is None:
                return
            if not keep_s.strip().lower().startswith("n"):
                self.mean_results.append((orientation, res))
                self._afficher(f"({len(self.mean_results)} mean result(s) kept in memory.)\n")

    def _format_mean_result_entry(self, index, orientation, res, code2=None):
        """Formate UNE entree de self.mean_results (ou lue depuis un
        .pmagani) pour l'affichage console - factorise entre "List
        results" et "Liste results file"."""
        site = mean_result_site(res.id)
        header = f"\n[{index}]"
        if code2:
            header += f" code2: {code2}"
        header += f" id: {res.id}"
        header += f"  (site: {site})" if site else "  (mixes more than one specimen/site prefix)"
        header += f"  orientation: {self._ORIENT_LABEL.get(orientation, '?')}"
        return header + "\n" + format_tsmean_box(res, orientation).rstrip("\n")

    def _current_pmagani_path(self):
        """Retourne self.ani_path s'il pointe vers un .pmagani (le seul
        format qui porte la section '#site mean tensor results') - None
        (+ avertissement) sinon. Utilise par "Select results"/"Liste
        results file"/"save results in File" - demande explicite
        utilisateur ("as the results are now in the same file, there is
        no need to open an other file when you select the results") :
        ces commandes operent desormais sur le fichier .pmagani DEJA
        ouvert (self.ani_path) plutot que d'en demander un autre."""
        if not self.ani_path:
            self._showwarning(
                "No file open",
                "Open a .pmagani file first (AMS Files > Open File .pmagani...).")
            return None
        if os.path.splitext(self.ani_path)[1].lower() != ".pmagani":
            self._showwarning(
                "Not a .pmagani file",
                f"The currently open file ({os.path.basename(self.ani_path)}) is a "
                "legacy .ANI file - site mean results need the new .pmagani format. "
                "Use \"Import legacy .ANI to .pmagani...\" first.")
            return None
        return self.ani_path

    def _save_mean_results_to_pmagani(self, results, source_label=""):
        """Ecrit chaque (orientation, TensorialMeanResult) de `results`
        dans la section '#site mean tensor results' du .pmagani
        COURAMMENT OUVERT (voir _current_pmagani_path et
        ams_selection.write_ani_mean_result) - demande explicite
        utilisateur ("as the results are now in the same file, there is
        no need to open an other file when you select the results") :
        specimen et mean partagent deja le meme fichier (voir le
        commentaire au-dessus de _ANI_SPECIMEN_HEADER cote
        STARpaleomag_Py/calcul), plus besoin d'en choisir un autre. Un
        resultat dont l'id ne se resout pas a un site unique (voir
        ams_selection.mean_result_site) ou qui est isotrope (pas d'axes
        propres, voir ams_selection._format_pmagani_mean_line) est ECARTE
        et journalise plutot qu'exporte sous un nom devine.

        `res.source_code2` ("RE"/"IM", voir ouvrir_tsmean_dialog qui
        separe deja Re/Im) FORCE le code2 ecrit pour CE resultat, sans
        redemander - demande explicite utilisateur ("force the code to
        Re and Im when saving the results in the file") : le prompt
        "Tensor type" n'est pose que s'il reste au moins un resultat SANS
        source_code2 connu (charge d'ailleurs, ou calcule par un chemin
        plus ancien), et ne s'applique alors qu'a ceux-la."""
        if not results:
            return
        path = self._current_pmagani_path()
        if path is None:
            return
        code2 = None
        if any(res.source_code2 is None for _o, res in results):
            code_s = self._console_input(
                "Tensor type for the remaining mean(s) without a known Re/Im source "
                "(A0=ATRM, F0=AARM, N0=AMS): ", "N0")
            if code_s is None:
                return
            code2 = code_s.strip().upper() or "N0"
        written, skipped = 0, []
        for orientation, res in results:
            site = mean_result_site(res.id)
            if site is None:
                skipped.append(f"{res.id or '(no id)'} (mixes more than one site)")
                continue
            try:
                write_ani_mean_result(
                    path, site, res.source_code2 or code2, res, orientation, info=source_label)
            except ValueError:
                skipped.append(f"{res.id} (isotropic)")
                continue
            written += 1
        msg = f"{written} mean tensor(s) saved to {path} (site mean section).\n"
        if skipped:
            shown = ", ".join(skipped[:10])
            more = "..." if len(skipped) > 10 else ""
            msg += f"Skipped {len(skipped)} result(s): {shown}{more}\n"
        self._afficher(msg)

    def ouvrir_selres_dialog(self):
        """Equivalent de `selres` (anisotropie.f:3050-3093) : charge des
        moyennes de site DEJA sauvegardees dans le .pmagani COURAMMENT
        OUVERT (voir _current_pmagani_path - demande explicite
        utilisateur "as the results are now in the same file, there is
        no need to open an other file when you select the results" :
        specimen et mean section vivent deja dans le meme fichier, plus
        besoin d'en ouvrir un autre), filtrees a l'orientation EXACTE de
        self.orientation, et laisse l'utilisateur choisir lesquelles
        ajouter a self.mean_results (liste de travail) - -1 = toutes,
        meme convention que le Fortran ("selection: (-1=all)").

        Filtre EXACT (iorient==iorient, meme regle que le Fortran
        d'origine) : IS et TC ne sont PAS regroupes - demande explicite
        utilisateur ("the select results still list both results in IS
        and TC during the selection. The TC and IS is not well
        constrained") suite a une version precedente qui les melangeait
        (rationale "ne differant que par une rotation de pendage connue"
        - abandonnee, pas assez solide pour garantir que les deux
        orientations restent comparables/interchangeables dans un
        calcul en aval)."""
        path = self._current_pmagani_path()
        if path is None:
            return
        all_entries = read_ani_mean_results_from_pmagani(path)
        if not all_entries:
            self._showwarning("No mean result", "No site mean tensor section found in this file.")
            return
        wanted_orientation = self.orientation.get()
        entries = [e for e in all_entries if e[0] == wanted_orientation]
        if not entries:
            self._showwarning(
                "No matching mean result",
                f"No mean result in {self._ORIENT_LABEL.get(self.orientation.get(), '?')} "
                f"coordinates in this file ({len(all_entries)} in other coordinate frame(s)).")
            return
        lines = [
            f"--- {len(entries)} mean tensor result(s) in "
            f"{self._ORIENT_LABEL.get(self.orientation.get(), '?')} coordinates "
            f"available in {os.path.basename(path)} "
            f"({len(all_entries) - len(entries)} in other coordinate frame(s) excluded) ---"]
        for i, (orientation, res, code2) in enumerate(entries, start=1):
            lines.append(
                f"[{i}] code2: {code2}  id: {res.id}  "
                f"orientation: {self._ORIENT_LABEL.get(orientation, '?')}")
        self._afficher("\n".join(lines) + "\n")
        sel_s = self._console_input("selection: comma-separated numbers, -1=all (0=cancel): ", "-1")
        if sel_s is None:
            return
        sel_s = sel_s.strip()
        if sel_s in ("0", ""):
            return
        if sel_s == "-1":
            chosen = list(range(1, len(entries) + 1))
        else:
            chosen = []
            for tok in sel_s.split(","):
                tok = tok.strip()
                if not tok:
                    continue
                try:
                    n = int(tok)
                except ValueError:
                    continue
                if 1 <= n <= len(entries):
                    chosen.append(n)
        for n in chosen:
            orientation, res, _code2 = entries[n - 1]
            self.mean_results.append((orientation, res))
        self._afficher(
            f"{len(chosen)} mean result(s) added to the in-memory list "
            f"({len(self.mean_results)} total).\n")

    def ouvrir_lisresmem_dialog(self):
        """Equivalent de `lisresmem`/"List results" (anisotropie.f:
        3094-3163) : liste self.mean_results (liste de travail en
        memoire), sans toucher a aucun fichier.

        Tableau compact (une ligne par site, voir
        ams_stats.format_lisresmem_table), PAS la boite detaillee
        (format_tsmean_box, celle affichee juste apres le calcul d'UN
        tenseur - toujours utilisee la, et par "Liste results file") -
        demande explicite utilisateur ("is it possible to use the
        original form of listing the mean tensor per site in memory
        (check the Fortran source) and not the detailed one page
        information") : verifie contre le Fortran d'origine (FORMAT 202),
        qui n'a jamais affiche de boite pour cette commande precise."""
        if not self.mean_results:
            self._showwarning("No mean result", "No mean tensor result in memory.")
            return
        header = f"--- {len(self.mean_results)} mean tensor result(s) in memory ---\n"
        self._afficher(header + format_lisresmem_table(self.mean_results))

    def ouvrir_initres_dialog(self):
        """Equivalent de `initres` (anisotropie.f:3244+) : remet
        self.mean_results a zero, apres confirmation (le Fortran demande
        aussi "mise a zero liste tenseurs moyens: Y/n")."""
        if not self.mean_results:
            self._afficher("Mean results list already empty.\n")
            return
        ans = self._console_input(
            f"Clear the {len(self.mean_results)} mean result(s) in memory? (y/N): ", "N")
        if ans is None or not ans.strip().lower().startswith("y"):
            return
        self.mean_results = []
        self._afficher("Mean results list cleared.\n")

    def ouvrir_delres_dialog(self):
        """Equivalent de `delres` (anisotropie.f:3165-3189) : supprime UNE
        entree de self.mean_results par son numero (voir "List results")."""
        if not self.mean_results:
            self._showwarning("No mean result", "No mean tensor result in memory.")
            return
        lines = [f"[{i}] {res.id}" for i, (_o, res) in enumerate(self.mean_results, start=1)]
        self._afficher("\n".join(lines) + "\n")
        idx_s = self._console_input("line to delete (0=cancel): ", "0")
        if idx_s is None:
            return
        try:
            idx = int(idx_s)
        except ValueError:
            self._showerror("Error", "Must be an integer.")
            return
        if idx <= 0 or idx > len(self.mean_results):
            return
        removed = self.mean_results.pop(idx - 1)
        self._afficher(f"Removed [{idx}]: {removed[1].id}\n")

    def ouvrir_lisresfich_dialog(self):
        """Equivalent de `lisresfich`/"Liste results file" (anisotropie.f:
        3190-3208) : liste, en LECTURE SEULE (ne modifie pas
        self.mean_results), la section mean deja sauvegardee du .pmagani
        COURAMMENT OUVERT (voir _current_pmagani_path - demande explicite
        utilisateur "as the results are now in the same file, there is
        no need to open an other file when you select the results").
        Pour charger certaines de ces entrees dans la liste de travail,
        voir "Select results"."""
        path = self._current_pmagani_path()
        if path is None:
            return
        entries = read_ani_mean_results_from_pmagani(path)
        if not entries:
            self._showwarning("No mean result", "No site mean tensor section found in this file.")
            return
        lines = [f"--- {len(entries)} mean tensor result(s) stored in {os.path.basename(path)} ---"]
        for i, (orientation, res, code2) in enumerate(entries, start=1):
            lines.append(self._format_mean_result_entry(i, orientation, res, code2=code2))
        self._afficher("\n".join(lines) + "\n")

    def ouvrir_sauveresfich_dialog(self):
        """Equivalent de `sauveresfich` (anisotropie.f:3210-3241) :
        sauvegarde self.mean_results (tenseurs moyens accumules en
        memoire - "Tensorial mean..." ou "Select results") dans la
        section mean du .pmagani COURAMMENT OUVERT (voir
        _current_pmagani_path) - demande explicite utilisateur ("on
        import, the results are not saved in the file pmagani" + "the
        menu within results are not implemented" + "as the results are
        now in the same file, there is no need to open an other file").
        A la difference du Fortran (deux listes distinctes
        amsres/amsresfichier), AMS_Py n'a qu'une seule liste de travail
        (self.mean_results) - c'est elle qui est sauvegardee en entier."""
        if not self.mean_results:
            self._showwarning("No mean result", "No mean tensor result in memory.")
            return
        self._save_mean_results_to_pmagani(self.mean_results, source_label="site mean (AMS_Py)")

    def ouvrir_corsondage_dialog(self):
        """Equivalent de `corsondage` (lect_asc.f:700-779) : correction
        d'orientation de carottier, deux methodes au choix comme le
        Fortran :
          1) "Orientation donnees de remanence avec foliation
             magnetique" : recalcule caz depuis la declinaison de l'axe
             k3 (mineur, normal au plan de foliation) du tenseur AMS de
             CHAQUE mesure de self.selection (voir
             ams_stats.principal_axes) - caz = 180 - dec(k3), ramene dans
             [0,360) (meme formule que le Fortran, ams(i).caz=180-
             axes(3,2,i)).
          2) "Orientation donnees d'anisotropie avec direction de
             remanence" : lit un fichier texte externe (2 colonnes : id,
             correction d'azimut en degres) et retranche cette correction
             du caz de chaque mesure de self.selection dont l'id
             correspond - comparaison de PREFIXE (simplification de la
             boucle caractere-par-caractere du Fortran, meme resultat
             pratique).
        Modifie self.selection EN PLACE (comme `ams(i).caz=...` cote
        Fortran) - pas de fichier de sortie separe."""
        if not self.selection:
            self._showwarning("No data", "No measurement selected.")
            return
        choice = self._console_input(
            "Orient remanence data with magnetic foliation (1) or "
            "orient AMS data with a known remanence direction (2)? : ", "1")
        if choice is None:
            return
        choice = choice.strip()
        if choice == "1":
            lines = ["--- Correction sondage (magnetic foliation) ---"]
            for m in self.selection:
                axes = principal_axes(m, self.orientation.get())
                _k3, dec3, _inc3 = axes[2]
                caz = 180.0 - dec3
                if caz > 360.0:
                    caz = 360.0 - caz
                elif caz < 0.0:
                    caz = 360.0 + caz
                lines.append(f"D  {m.id}  cin={m.cin:.1f}  caz={caz:.1f} (was {m.caz:.1f})")
                m.code2 = "C0"
                m.etape = 0
                m.caz = caz
            self._afficher("\n".join(lines) + "\n")
        elif choice == "2":
            path = filedialog.askopenfilename(
                title="Azimuth correction file (id, azimuth correction)",
                filetypes=[("Text", "*.txt *.asc"), ("All files", "*.*")])
            if not path:
                return
            corrections = []
            with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                    try:
                        corrections.append((parts[0], float(parts[1])))
                    except ValueError:
                        continue
            lines = ["--- Correction sondage (known remanence direction) ---"]
            applied = 0
            for m in self.selection:
                for echan, carcaz in corrections:
                    if echan[:len(m.id)] == m.id or m.id[:len(echan)] == echan:
                        before = m.caz
                        m.caz = m.caz - carcaz
                        lines.append(f"{m.id}  azimuth before: {before:.1f}  after: {m.caz:.1f}")
                        applied += 1
                        break
            lines.append(f"\n({applied}/{len(self.selection)} measurement(s) corrected.)")
            self._afficher("\n".join(lines) + "\n")
        else:
            self._showerror("Error", "Choice must be 1 or 2.")

    def ouvrir_mds_dialog(self):
        """Equivalent de `mds` (anisotropie.f:3560-3623) : moyenne
        arithmetique et geometrique (+ ecart-type geometrique) de
        `axes(2,1,i)*1e-5` (l'axe k2 du tenseur DE-NORMALISE, PAS `s`
        directement - `editams` reconstruit le tenseur physique via
        `k_ij*|s|` avant `elprop`, donc axes(2,1,i) est deja a l'echelle
        physique de `s`, d'ou le facteur 1e-5 - transcrit tel quel, meme si
        utiliser k2 plutot que la susceptibilite bulk (trace/3) comme
        "susceptibilite moyenne" est une convention propre a cette
        routine). Necessite au moins 3 mesures."""
        if len(self.selection) < 3:
            self._showwarning("Not enough data", "Mean susceptibility needs at least 3 measurements.")
            return
        import numpy as np
        suscint = []
        for m in self.selection:
            a = np.array([
                [m.k11 * abs(m.s), m.k12 * abs(m.s), m.k13 * abs(m.s)],
                [m.k12 * abs(m.s), m.k22 * abs(m.s), m.k23 * abs(m.s)],
                [m.k13 * abs(m.s), m.k23 * abs(m.s), m.k33 * abs(m.s)],
            ])
            av = np.linalg.eigvalsh(a)[::-1]  # decroissant
            suscint.append(av[1] * 1.0e-5)  # k2
        res = mean_susceptibility(suscint)
        if res is None:
            self._showwarning("Not enough data", "Mean susceptibility needs at least 3 measurements.")
            return
        self._afficher(
            f"site: {self.selection[-1].id}  susceptibility  Nb: {res['n']:4d}  "
            f"Geom mean.: {res['geom_mean']:.3E}  Arith mean: {res['arith_mean']:.3E}\n"
        )

    def ouvrir_bootstrap_dialog(self):
        """Equivalent de `bootams` (menu Calcul > ellipses Bootstrap) - le
        Fortran d'origine est un STUB VIDE (`subroutine bootams; return;
        end`, AMS_OSX_x.f95:391-393), il n'y a donc rien a porter depuis le
        source. A la demande de l'utilisateur, implementation neuve basee
        sur `pmagpy.pmag.s_boot`/`sbootpars` (paquet PyPI installe -
        `pip install pmagpy`) : bootstrap non parametrique par defaut sur
        la liste de mesures couramment en memoire (`self.selection`), avec
        ajustement d'une distribution de Kent (zeta/eta) sur les vecteurs
        propres rechantillonnes - methode differente de l'ellipse de
        Jelinek deja portee (`Calcul > Tensorial mean`), voir
        ams_bootstrap.py. Necessite au moins 3 mesures."""
        if len(self.selection) < 3:
            self._showwarning(
                "Not enough data", "Bootstrap mean needs at least 3 measurements.")
            return
        nb_s = self._console_input("number of bootstrap draws : ", "1000")
        if nb_s is None:
            return
        try:
            nb = int(nb_s)
        except ValueError:
            self._showerror("Error", "Must be an integer.")
            return
        par_s = self._console_input("parametric bootstrap ? (y/N) : ", "N")
        if par_s is None:
            return
        parametric = par_s.strip().lower().startswith("y")

        res = compute_bootstrap_mean(
            self.selection, self.orientation.get(), nb=nb, parametric=parametric)
        if res is None:
            self._showwarning(
                "Not enough data", "Bootstrap mean needs at least 3 measurements.")
            return
        self._afficher(format_bootstrap_result(res))
        self._last_bootstrap = res

        plot_s = self._console_input(
            "Plot bootstrap ? cloud(0) ellipse(1) both(2) none(n) : ", "2")
        if plot_s is None:
            return
        plot_s = plot_s.strip().lower()
        if plot_s.startswith("n"):
            return
        self._bootstrap_show_cloud = plot_s in ("0", "2")
        self._bootstrap_show_ellipse = plot_s in ("1", "2")
        self._current_graphic = "bootstrap"
        self._refresh_current_graphic()


def main():
    root = tk.Tk()
    # Force l'encodage systeme de Tcl a utf-8 - demande explicite
    # utilisateur ("l'appli installee ne fonctionne pas tres bien, par
    # exemple probleme de texte Latin lors de l'importation. Pas de pb
    # depuis le terminal") : Tcl/Tk devine son "system encoding" depuis
    # LANG/LC_ALL au demarrage - un Terminal herite la locale du shell
    # (LANG deja UTF-8), un .app lance depuis le Finder/Dock n'a
    # generalement AUCUNE locale definie, Tcl se rabat alors sur un
    # encodage non-UTF-8 pour les widgets texte, meme si la chaine
    # Python elle-meme est deja correctement decodee (accents mal
    # affiches dans l'UI, pas une erreur de lecture de fichier).
    root.tk.call("encoding", "system", "utf-8")
    AmsApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
