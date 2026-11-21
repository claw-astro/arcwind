---
name: Bug report
about: Something does not work, or gives a wrong result
labels: bug
---

**What happened**
A clear description of the problem, with the full error message if there is one.

**What you expected**

**How to reproduce**
The configuration (TOML) and the Python lines or `arcwind` command you ran.
If possible, reproduce it on a synthetic cube (`arcwind.mock.MockWind.elias20_like()`).
Please don't attach large data cubes; describe the cube (size, channel width, beam) instead.

**Versions**
- arcwind: (`python -c "import arcwind; print(arcwind.__version__)"`)
- Python:
- OS:
