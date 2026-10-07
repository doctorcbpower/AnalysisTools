#!/usr/bin/env python3
"""Diagnostic: FOF membership, R200c crossings and pericentres of selected subhaloes; FOF-envelope statistics.

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/10_diag_group_membership.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import logging, numpy as np
from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.orbits import *
from analysistools.orbits.kinematics import host_reference
logging.basicConfig(level=logging.ERROR)
mt=MergerTreeTools(TREE,'SubFind',comoving=True,little_h=True,native_includes_h=True,native_is_comoving=True,loglevel=logging.ERROR)
cosmo=FlatLCDM.from_tree(mt); d=mt.data; lut=d.lookup
units=TrackUnits(comoving=True,little_h=True,h=0.678,length_kpc=1000.0,mass_msun=1e10,velocity='physical')
host=mt.get_track(0,122); members=np.array(mt.get_group_members(0,122))
m0=np.array([d.SubhaloMass[lut[(122,int(m))]] for m in members]); top=members[np.argsort(m0)[::-1][:10]]
cfg=OrbitConfig(radius_key='R200c_from_M200')
tr=OrbitTracker(cosmo,units,cfg,boxsize=mt.BoxSize)
allt=[mt.get_track(int(m),122) for m in members]
cat=tr.run(allt,host,ids=members.tolist())
rref=host_reference(host,units,'R200c_from_M200',1.0,None,cosmo)
print('R200c(z=0) = %.1f kpc'%rref[-1])
pos={int(m):i for i,m in enumerate(members)}
for sid in (3,4,7,9):
    i=pos[sid]; r=cat.record(i); S=r.series; sc=r.scalars; t=allt[i]
    ratio=S['r']/S['r_ref']
    # FOF membership history from the tree (IsSubhalo = not first in its FOF group)
    fof=t.IsSubhalo
    ins=[e for e in r.events if e.kind==CROSS_IN]; outs=[e for e in r.events if e.kind==CROSS_OUT]
    print(f"\nSubhaloNr {sid}: r/R200c at z=0 = {ratio[-1]:.2f}; min r/R200c over history = {ratio.min():.2f} at a={S['a'][ratio.argmin()]:.3f} (r={S['r'][ratio.argmin()]:.0f} kpc)")
    print(f"   FOF satellite at z=0: {bool(fof[-1])};  first FOF-satellite snapshot a={t.Time[np.argmax(fof)]:.3f}; fraction of snapshots as FOF sat = {fof.mean():.2f}")
    print(f"   R200c crossings in at a={[round(e.a,3) for e in ins]}  out at a={[round(e.a,3) for e in outs]}  -> backsplash={bool(sc['is_backsplash'])}, inside_now={bool(sc['inside_now'])}")
    print(f"   pericentres (kpc, a): {[(round(e.r,0),round(e.a,3)) for e in r.events if e.kind==PERICENTRE]}")
# FOF envelope at z=0: r/R200c for all members
rf=np.array([cat.record(i).series['r'][-1] for i in range(len(members))])/rref[-1]
print('\nAll %d z=0 FOF members: r/R200c percentiles 50/90/99/max ='%len(members),np.percentile(rf,[50,90,99,100]).round(2),'; fraction outside R200c = %.2f'%np.mean(rf>1))
big=np.array([units.mass(d.SubhaloMass[lut[(122,int(m))]]) for m in members])>1e9
print('For the %d members > 1e9 Msun: fraction outside R200c = %.2f, max r/R200c = %.2f'%(big.sum(),np.mean(rf[big]>1),rf[big].max()))
