## AnalysisTools - scripting tools for analysing astrophysical simulations

These assume simulations are GADGET2/3/4, Arepo, or SWIFT-based. It is straightforward(-ish) to add additional functionality for other codes.

Install using `pip install -e .`.

Load as python module, `import analysistools`

You can then load data using,
```
snap = analysistools.SnapshotTools(snapfileformat="HDF5");
data = snap.read_snapshot("./data/snap_0031", convention="SWIFT")  # Omit the suffix of the snapshot
```
and then access the data with `snap.pos`, `snap.vel`, etc...

You can separate particles by type,
```
snap.LoadParticlesByType(part_type="all")    # Can be all, gas, star, bh
```
and then access with `snap.dm.pos`, `snap.dm.vel`.

You can apply conversions to comoving/physical and to length/h or simply length using
```
snap.UnitConversion('convert_to_comoving','convert_to_per_littleh')    # Alterantives are convert_to_physical and convert_to_littleh
```

If you want to write a snapshot (say, converting between formats, or only a subset of particles), you can use,
```
idx_type=snap.ParticleOffsetsByType(snap.num_part_total)
import numpy as np
snap.write_snapshot("./data/snapshot_031",idx=np.arange(snap.pids.size),idx_type=idx_type,convention="arepo",
                    blocks_to_write=['pos','vel'])
```


