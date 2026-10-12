# AMS_Py

Analysis of anisotropy of magnetic susceptibility (AMS) data: import of AGICO kappabridge files, stereoplots, tensor means and statistics, selection for MagIC export. Works with the `.prmag` / `.pmagani` files of [STARpaleomag_Py](https://github.com/pierrickr35-breton/STARpaleomag_Py).

AMS_Py is a plain **Python 3** program (Tkinter + Matplotlib). It runs on **Windows, Linux and macOS**. For Macs with an Apple Silicon processor, a ready-to-use application is also provided.

## Running on any platform (Windows, Linux, macOS)

1. Install [Python 3](https://www.python.org/downloads/) if it is not already on your computer. On Linux, Tkinter may be a separate package (for example `sudo apt install python3-tk`).
2. Get the program: **[download the source code (zip)](https://github.com/pierrickr35-breton/AMS_Py/archive/refs/heads/main.zip)** and unzip it, or clone this repository with git.
3. Install the required libraries (once):

   ```bash
   pip install numpy matplotlib pmagpy
   ```

4. Start the program from its folder:

   ```bash
   python3 app.py
   ```

   (on Windows: `python app.py`)

To update later, download the source code again (or `git pull`): the libraries do not need to be reinstalled.

## macOS: ready-to-use application

For Macs with an **Apple Silicon** processor (M1 and later), a standalone application is available: Python and all libraries are included, there is nothing else to install.

1. **[Go to the latest release](https://github.com/pierrickr35-breton/AMS_Py/releases/latest)** and, under *Assets*, click `AMS_Py.zip`.
2. Double-click the zip to unpack it, and move `AMS_Py.app` to your *Applications* folder.
3. The application is not signed with an Apple developer certificate, so macOS blocks it the first time. Try to open it once, then go to *System Settings > Privacy & Security*, scroll down and click **Open Anyway**. This is needed only once.

   If macOS says the application "is damaged and can't be opened", open the Terminal and type:

   ```bash
   xattr -dr com.apple.quarantine /Applications/AMS_Py.app
   ```

On an Intel Mac, use the Python instructions above.

## User Guide

In the program: *Help > User Guide*. The same guide is in this repository: `help/AMS_Py_Guide.html`.
