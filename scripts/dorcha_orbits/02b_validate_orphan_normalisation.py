#!/usr/bin/env python3
"""Same truncation test comparing the NFW normalisations 'm200' and 'vmax'.

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/02b_validate_orphan_normalisation.py
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

from analysistools.orbits.potential import NFWPotential, _nfw_m
from analysistools.orbits.cosmology import G_KPC_KMS2_MSUN as G
from analysistools.orbits.kinematics import host_reference
a=host.Time; rmax=units.physical_kpc(host.extra['VmaxRad'],a); vmax=host.extra['Vmax']; r200=host_reference(host,units,'R200c_from_M200',1.0,None,cosmo)
def vmax_norm_pot():
    ts,ps=[],[]
    for i in range(len(a)):
        if rmax[i]<=0 or vmax[i]<=0: continue
        rs=rmax[i]/2.1626; m_rmax=vmax[i]**2*rmax[i]/G; norm=m_rmax/float(_nfw_m(2.1626))
        c=r200[i]/rs; m200=norm*float(_nfw_m(c)); ps.append(NFWPotential(m200,r200[i],rs)); ts.append(float(cosmo.age(a[i])))
    return TimeDependentPotential(ts,ps,r_grid=np.geomspace(0.5,3e3,160))
pots={'M200-normalised (current)':pot,'Vmax-normalised':vmax_norm_pot()}
cfg=OrbitConfig(radius_key='R200c_from_M200',orphan_max_r_over_ref=np.inf)
truth=OrbitTracker(cosmo,units,cfg,potential=pot,boxsize=mt.BoxSize).run(tracks,host,ids=members); rtrue=truth.scalar('r_final')
for cut in (100,80):
    tt=[];keep=[]
    for i,t in enumerate(tracks):
        if t.SnapNum[0]<=cut-5 and t.SnapNum[-1]>=122: tt.append(trunc(t,cut)); keep.append(i)
    keep=np.array(keep)
    for name,pp in pots.items():
        pred=OrbitTracker(cosmo,units,cfg,potential=pp,boxsize=mt.BoxSize).run(tt,host,ids=[members[i] for i in keep])
        rp=pred.scalar('r_final'); lastok=[np.flatnonzero(~pred.record(i).series['extrapolated'])[-1] for i in range(len(tt))]
        r0=np.array([pred.record(i).series['r'][k] for i,k in enumerate(lastok)]); rref=np.array([pred.record(i).series['r_ref'][k] for i,k in enumerate(lastok)])
        sel=r0<rref; e=np.log10(rp[sel]/rtrue[keep][sel])
        print(f'cut snap {cut}  {name:28s} N={sel.sum():4d} med dlog10r={np.median(e):+.3f}  16-84%=[{np.percentile(e,16):+.2f},{np.percentile(e,84):+.2f}]  within30%={np.mean(abs(rp[sel]/rtrue[keep][sel]-1)<0.3):.2f}  pred<R200 {np.mean(rp[sel]<273):.2f} (true {np.mean(rtrue[keep][sel]<273):.2f})')
