#!/usr/bin/env python3
"""Continue lost (branch-ending) subhaloes of the host as orphans and report where they end up; shows why continuation is gated by orphan_max_r_over_ref.

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/03_orphan_fringe_check.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import time, logging, numpy as np
from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.orbits import *
logging.getLogger('analysistools.orbits').setLevel(logging.INFO)
logging.basicConfig(level=logging.WARNING)
mt=MergerTreeTools(TREE,'SubFind',comoving=True,little_h=True,native_includes_h=True,native_is_comoving=True,loglevel=logging.ERROR)
d=mt.data; cosmo=FlatLCDM.from_tree(mt)
units=TrackUnits(comoving=True,little_h=True,h=0.678,length_kpc=1000.0,mass_msun=1e10,velocity='physical')
host=mt.get_track(0,122)
pot=tree_nfw_potential(host,cosmo,units,radius_key='R200c_from_M200')
H=__import__('h5py').File(TREE,'r')['TreeHalos']
desc=H['TreeDescendant'][()]; snap=np.asarray(d.SnapNum); grp=np.asarray(d.GrpNr)
host_grp=host.extra['GroupNr']; hs=host.SnapNum
# nodes in the host's FOF group at their snapshot, branch ends (no descendant), before z=0
sel=[]
for s,g in zip(hs,host_grp):
    if s>=122 or s<60: continue
    idx=np.flatnonzero((snap==s)&(grp==g)&(desc==-1)); sel.extend(idx.tolist())
sel=np.array(sel); print('lost (branch-ending) subhaloes in host FOF, snaps 60-121:',sel.size)
rng=np.random.default_rng(0); sel=rng.choice(sel,size=min(1500,sel.size),replace=False)
tracks=[mt.get_track(int(d.SubhaloNr[i]),int(snap[i])) for i in sel]
for df in (False,):
    tr=OrbitTracker(cosmo,units,OrbitConfig(radius_key='R200c_from_M200',dynamical_friction=df),potential=pot,boxsize=mt.BoxSize)
    t0=time.time(); cat=tr.run(tracks,host,ids=sel.tolist()); dt=time.time()-t0
    S={k:cat.scalar(k) for k in SCALARS if k!='L_infall'}
    last=np.array([t.SnapNum[-1] for t in tracks]); m_last=np.array([units.mass(t.Mass[-1]) for t in tracks])
    # r at last resolved snapshot
    r_last=np.array([cat.record(i).series['r'][np.flatnonzero(~cat.record(i).series['extrapolated'])[-1]] for i in range(len(tracks))])
    print(f'{dt:.1f}s  tree-ended={S["is_orphan_final"].sum()}  continued={S["orphan_continued"].sum()}/{len(tracks)}')
    c=S['orphan_continued']
    print('  CONTINUED subset: r_last',np.percentile(r_last[c],[5,50,95]).round(1),'r_final',np.percentile(S['r_final'][c],[5,50,95]).round(1),' frac r_final<R200c',np.mean(S['r_final'][c]<273).round(3),' <20kpc',np.mean(S['r_final'][c]<20).round(3))
    print('  r_last_resolved [kpc] 5/50/95:',np.percentile(r_last,[5,50,95]).round(1),'  r_final:',np.percentile(S['r_final'],[5,50,95]).round(1))
    print('  r_final<R200c(z0=273):',np.mean(S['r_final']<273).round(3),' r_final<20kpc:',np.mean(S['r_final']<20).round(3), ' non-finite:',(~np.isfinite(S['r_final'])).sum())
    print('  last resolved snap distribution (pct):',np.percentile(last,[5,50,95]), ' mass_last (Msun) median',np.median(m_last))
