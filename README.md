NHC_project
===========

Small molecular dynamics (MD) code using a Nose–Hoover chain (NHC) thermostat.

Repository structure
--------------------
- `src/` : C++ source files
- `include/` : project headers
- `Makefile` : build script (produces the `md` executable)
- `data.py` : Python script to load `output.log` and produce plots
- `input-16.in`, `input-500.in` : example input parameter files
- `coords16.data`, `coords500.data` : example coordinate files

Quick start (Linux / WSL / macOS)
--------------------------------
Requirements:
- g++ (supporting C++11)
- make
- Python 3 with `numpy` and `matplotlib`

Build:

```bash
make
```

Run:

```bash
./md input-16.in
```

The program writes `output.log` and calls `data.py` to generate `md_energy_temp.png`.

Plotting manually:

```bash
python3 data.py
```

Windows
-------
You can build and run on Windows using one of these approaches:

1. WSL (recommended): install Ubuntu/WSL, then follow the Linux instructions above.

2. MSYS2 / MinGW: install MSYS2, use `pacman` to install mingw-w64 toolchain and Python; build with `mingw32-make` or `make`.

3. Visual Studio: create a project or CMake target and add `src/` files; build with the IDE.

Python dependencies
-------------------
Install Python 3 and then:

```bash
pip3 install numpy matplotlib
```

Notes
-----
- Consider adding `.gitignore` to exclude binary artifacts (`bin/`, `*.o`, `md`, `output.log`, `*.png`).
- This repo was recently reset to keep the current working tree; remote history may have been replaced.

If you want, I can add a `.gitignore` and a Chinese README (`README-zh.md`) and push them as well.

