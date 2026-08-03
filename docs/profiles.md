# Profiles and Mass Functions

`ProfileTools` computes radial/vertical density and kinematic profiles; `MassFunctionTools` (same module) computes halo/galaxy mass functions.

```python
from analysistools.profile_tools import ProfileTools

pt = ProfileTools(numbins=40)

rho = pt.volume_density(pos, mass, centre, rmin=0.1, rmax=300)
sigma = pt.surface_density(pos, mass, centre, rmin=0.1, rmax=50)
hz = pt.scale_height(pos, mass, centre, rmin=0.1, rmax=50)

pt.plot(rho, "density", ylabel=r"$\rho(r)$")
```

For the dataset-aware equivalents (`plt2.profile(...)`, `plt2.mass_function(...)`) that work directly on `Dataset` objects, see [unified_interface.md](unified_interface.md).
