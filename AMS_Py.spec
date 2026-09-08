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

a = Analysis(
    ['app.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=[],
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
)
