# AMS_Py

Analysis of anisotropy of magnetic susceptibility (AMS) data: import of AGICO kappabridge files, stereoplots, tensor means and statistics, selection for MagIC export. Works with the `.prmag` / `.pmagani` files of [STARpaleomag_Py](https://github.com/pierrickr35-breton/STARpaleomag_Py).

## Download (macOS standalone application)

**[Download the latest version](https://github.com/pierrickr35-breton/AMS_Py/releases/latest)** — on the release page, under *Assets*, click `AMS_Py.zip`.

- This is a **standalone macOS application**: Python and all libraries are included, there is nothing else to install.
- It runs on Macs with an **Apple Silicon** processor (M1 and later) only. On an Intel Mac, Windows or Linux, run the program from the source code (see below).
- Do **not** use the "Source code (zip)" / "Source code (tar.gz)" links of the release page unless you want the Python source: they do not contain the application.

### Installation

1. Double-click `AMS_Py.zip` to unzip it, and move `AMS_Py.app` to your *Applications* folder.
2. The application is not signed with an Apple developer certificate, so macOS blocks it the first time. Try to open it once, then go to *System Settings > Privacy & Security*, scroll down and click **Open Anyway**. This is needed only once.

   If macOS says the application "is damaged and can't be opened", open the Terminal and type:

   ```bash
   xattr -dr com.apple.quarantine /Applications/AMS_Py.app
   ```

## User Guide

In the application: *Help > User Guide*. The same guide is in this repository: `help/AMS_Py_Guide.html`.

## Running from the source code (any platform)

Requires Python 3 with Tkinter, and: numpy matplotlib pmagpy.

```bash
pip install numpy matplotlib pmagpy
python3 app.py
```
