# FDM (Fuzzy Dark Matter) Field Snapshots

`FDMFieldTools` is a mesh-field analogue to `SnapshotTools`, for reading FDM wavefunction field snapshots (not particle data):

```python
from analysistools.fdm_field_tools import FDMFieldTools

field = FDMFieldTools("fdm_field_010.hdf5")

field.psi       # complex128 array, shape (N, N, N)
field.density   # |psi|^2, real array, shape (N, N, N)
field.N         # grid size
field.L         # box size (code length units)
field.mass      # mc^2, eV
```
