# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_data_files

# ams_bootstrap.py imports pmagpy.pmag - meme raison que STARpaleomag_Py.spec :
# pmagpy expedie des fichiers de donnees hors .py (field_models/,
# data_model/) que PyInstaller ne detecte pas automatiquement.
datas = collect_data_files('pmagpy')
# Guide utilisateur statique (Help > User Guide, voir app._resource_path/
# ouvrir_user_guide) - meme raison que STARpaleomag_Py.spec : doit etre
# EXTRAIT sous le meme nom de dossier ('help/') pour que _resource_path
# (sys._MEIPASS + 'help' + nom de fichier) le retrouve une fois empaquete.
datas += [('help', 'help')]

# matplotlib charge le backend SVG dynamiquement au moment de
# fig.savefig(path, format="svg") (Export SVG...), PAS via un `import`
# statique visible dans le source - PyInstaller ne le detecte donc pas
# tout seul et l'app packagee echoue a l'export avec "No module named
# matplotlib.backends.backend_svg" (repere sur un vrai .app construit,
# fonctionnait depuis les sources). backend_pdf ajoute par precaution/
# coherence avec les autres .spec de ce projet, meme si non utilise ici.
hiddenimports = ['matplotlib.backends.backend_svg', 'matplotlib.backends.backend_pdf']

a = Analysis(
    ['app.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='AMS_Py',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='AMS_Py',
)
app = BUNDLE(
    coll,
    name='AMS_Py.app',
    icon=None,
    bundle_identifier=None,
    # LANG/LC_ALL explicites - un .app lance depuis le Finder/Dock n'a
    # generalement AUCUNE locale definie (contrairement a un Terminal,
    # qui herite celle du shell) - demande explicite utilisateur
    # ("l'appli installee ne fonctionne pas tres bien, par exemple
    # probleme de texte Latin lors de l'importation. Pas de pb depuis le
    # terminal"). Complement de root.tk.call("encoding", "system",
    # "utf-8") dans app.py (celui-ci force Tcl/Tk specifiquement ; ceci
    # couvre tout le reste du processus - locale.*, etc.).
    info_plist={'LSEnvironment': {'LANG': 'en_US.UTF-8', 'LC_ALL': 'en_US.UTF-8'}},
)
