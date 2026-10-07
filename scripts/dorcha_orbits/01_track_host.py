#!/usr/bin/env python3
"""Track every z=0 FOF member of the host (group 0) with the orbit tool and print summary statistics; writes the sidecar orbits_host0.h5.

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/01_track_host.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import time, logging, numpy as np
logging.basicConfig(level=logging.INFO)
from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.orbits import *
t0=time.time()
mt=MergerTreeTools(TREE,'SubFind',comoving=True,little_h=True,
                   native_includes_h=True,native_is_comoving=True)
print('loaded',round(time.time()-t0,1),'s; BoxSize',mt.BoxSize)
d=mt.data
snap=np.asarray(d.SnapNum); first=np.asarray(d.TreeFirstHaloInFOFgroup)
cen=np.flatnonzero((snap==122)&(first==np.arange(snap.size)))
host_idx=cen[np.argmax(np.asarray(d.SubhaloMass)[cen])]
host_sub=int(d.SubhaloNr[host_idx]); grp=int(d.GrpNr[host_idx])
print('host SubhaloNr',host_sub,'GroupNr',grp,'M200 (1e10Msun/h)',d.GrpM200[host_idx])
host=mt.get_track(host_sub,122)
members=mt.get_group_members(grp,122)
print('satellites in FOF group at z=0:',len(members))
tracks=[mt.get_track(int(m),122) for m in members]
cosmo=FlatLCDM.from_tree(mt)
units=TrackUnits(comoving=True,little_h=True,h=0.678,length_kpc=1000.0,mass_msun=1e10,velocity='physical')
cfg=OrbitConfig(radius_key='R200c_from_M200')
pot=tree_nfw_potential(host,cosmo,units,radius_key='R200c_from_M200')
print(pot.t.size,'NFW epochs; first/last:'); 
tr=OrbitTracker(cosmo,units,cfg,potential=pot,boxsize=mt.BoxSize)
t0=time.time(); cat=tr.run(tracks,host,ids=[int(m) for m in members]); print('tracked in',round(time.time()-t0,1),'s')
cat.write_hdf5(os.path.join(OUT_DIR,'orbits_host0.h5'))
import pickle
S={k:cat.scalar(k) for k in SCALARS if k!='L_infall'}
r200=r_ = None
from analysistools.orbits.kinematics import host_reference
rref=host_reference(host,units,'R200c_from_M200',1.0,None,cosmo)
print('host R200c z=0 [kpc]:',rref[-1],' z=1:',rref[np.argmin(abs(host.Redshift-1))])
print('host Vmax/Rmax z0:',host.extra['Vmax'][-1],host.extra['VmaxRad'][-1]*1000,'kpc(comoving=phys at z0)')
for k in ['born_inside','n_infalls','is_backsplash','inside_now','preprocessed','n_pericentres','is_orphan_final','n_snapshot_gaps']:
    print(k, int(S[k].sum()) if S[k].dtype==bool else np.bincount(S[k].astype(int))[:8], )
fin=S['r_final']; print('r_final: median',np.median(fin),'max',fin.max(),' frac<R200c',np.mean(fin<rref[-1]))
for k in ['r_peri_min','r_apo_last','ecc_last','z_infall','mass_infall','period_peri','r_peri_EL','r_apo_EL','period_EL','E_infall']:
    v=S[k]; g=np.isfinite(v); print(f'{k:12s} n={g.sum():4d}', (np.percentile(v[g],[5,50,95]).round(3) if g.any() else ''))
