#!/usr/bin/env python3
"""Truncation test of orphan continuation: cut satellites the tree follows to z=0 at several epochs, continue them as orphans, compare with the tree's own z=0 radius (M200-normalised NFW).

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/02a_validate_orphans_truncation.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import logging, copy, numpy as np
from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.merger_tree_types import HaloTrack
from analysistools.orbits import *
logging.basicConfig(level=logging.ERROR)
mt=MergerTreeTools(TREE,'SubFind',comoving=True,little_h=True,native_includes_h=True,native_is_comoving=True,loglevel=logging.ERROR)
cosmo=FlatLCDM.from_tree(mt)
units=TrackUnits(comoving=True,little_h=True,h=0.678,length_kpc=1000.0,mass_msun=1e10,velocity='physical')
host=mt.get_track(0,122)
pot=tree_nfw_potential(host,cosmo,units,radius_key='R200c_from_M200')
members=mt.get_group_members(0,122)
tracks=[mt.get_track(int(m),122) for m in members]
def trunc(t,snap):
    k=t.SnapNum<=snap
    return HaloTrack(halo_id=t.halo_id,query_snapnum=int(snap),treefileformat=t.treefileformat,SnapNum=t.SnapNum[k],Redshift=t.Redshift[k],Time=t.Time[k],Mass=t.Mass[k],Pos=t.Pos[k],Vel=t.Vel[k],IsSubhalo=t.IsSubhalo[k],HostID=t.HostID[k],extra={kk:v[k] for kk,v in t.extra.items()})
cfg=OrbitConfig(radius_key='R200c_from_M200',orphan_max_r_over_ref=np.inf)
tr=OrbitTracker(cosmo,units,cfg,potential=pot,boxsize=mt.BoxSize)
truth=tr.run(tracks,host,ids=members)          # full tree orbits (reference)
rtrue=truth.scalar('r_final')
for cut,label in [(110,'cut z~0.1 (~1.3 Gyr before end)'),(100,'cut z~0.3'),(80,'cut z~0.8'),(65,'cut z~1.5')]:
    tt=[]; keep=[]
    for i,t in enumerate(tracks):
        if t.SnapNum[0]<=cut-5 and t.SnapNum[-1]>=122:
            k=t.SnapNum<=cut
            # inside R200c at the cut
            tt.append(trunc(t,cut)); keep.append(i)
    keep=np.array(keep)
    pred=tr.run(tt,host,ids=[members[i] for i in keep])
    rp=pred.scalar('r_final'); r0=np.array([pred.record(i).series['r'][np.flatnonzero(~pred.record(i).series['extrapolated'])[-1]] for i in range(len(tt))])
    rr=np.array([pred.record(i).scalars['r_final'] for i in range(len(tt))])
    rt=rtrue[keep]; rcut_ref=np.array([pred.record(i).series['r_ref'][np.flatnonzero(~pred.record(i).series['extrapolated'])[-1]] for i in range(len(tt))])
    for sel_name,sel in [('inside R200c at cut',r0<rcut_ref),('r<0.5 R200c at cut',r0<0.5*rcut_ref)]:
        e=np.log10(rp[sel]/rt[sel]); 
        print(f'{label:34s} {sel_name:20s} N={sel.sum():4d}  log10(r_pred/r_true) med={np.median(e):+.3f} 16-84%=[{np.percentile(e,16):+.2f},{np.percentile(e,84):+.2f}]  frac within 30% = {np.mean(abs(rp[sel]/rt[sel]-1)<0.3):.2f}   true<R200:{np.mean(rt[sel]<273):.2f} pred<R200:{np.mean(rp[sel]<273):.2f}')
