## AnalysisTools - scripting tools for analysing astrophysical simulations

These assume simulations are GADGET2/3/4, Arepo, or SWIFT-based. It is straightforward(-ish) to add additional functionality for other codes.

Install using `pip install -e .`.

Load as python module, `import analysistools`

You can then load data using,
```
snap = analysistools.snapshot_tools.SnapshotTools(snapfileformat="HDF5")
data = snap.read_snapshot("./snapshot_122", convention="GADGET4") # Omit the suffix of the snapshot
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
idx = np.where(mask)[0]
ptype = np.ones(len(idx), dtype=np.int64)

snap.write_snapshot(filename='/Users/00075868/CurrentWork/Dorcha/zoom/snap_122.cube',
                    file_format='HDF5',
                    convention='GADGET4',
                    idx=idx,
                    idx_type=ptype)
```


