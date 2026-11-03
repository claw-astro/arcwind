# Schematic figures of the Elias 2-20 paper

These scripts draw the model-definition schematics (`make_model_definitions.py`,
`make_model_definitions_v2.py`) and the 3-D ring schematic (`make_ring_schematic.py`).
They read the saved fit `<output.dir>/broken-power-law.json`; `wkcompat.py` adapts them
to diskwind.

```bash
python make_model_definitions.py --config ../configs/elias20.toml
python make_model_definitions_v2.py --config ../configs/elias20.toml
python make_ring_schematic.py --config ../configs/elias20.toml
```
