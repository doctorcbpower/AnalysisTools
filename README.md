## AnalysisTools

Install using `pip install -e .`.

Load as python module, `import analysistools`

You can then load data using,
```
snap=analysistools.SnapshotTools(snapfilename='data/snap_0031',snapfileformat='HDF5',convention='SWIFT')
snap.ReadSnapshot()
```


