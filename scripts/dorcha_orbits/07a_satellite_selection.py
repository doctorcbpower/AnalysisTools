#!/usr/bin/env python3
"""Satellite selection inside R200c at z=0: counts, FOF vs other-FOF members.

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/07a_satellite_selection.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import logging, numpy as np
from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.orbits import *
logging.basicConfig(level=logging.ERROR)
mt=MergerTreeTools(TREE,'SubFind',comoving=True,little_h=True,native_includes_h=True,native_is_comoving=True,loglevel=logging.ERROR)
cosmo=FlatLCDM.from_tree(mt); d=mt.data
units=TrackUnits(comoving=True,little_h=True,h=0.678,length_kpc=1000.0,mass_msun=1e10,velocity='physical')
host=mt.get_track(0,122)
sel=select_satellites_subfind(mt,host,cosmo,units,122)
print('R200c = %.1f kpc ; subhaloes inside R200c (excl. host central): %d'%(sel.r_ref,len(sel)))
print('  in host FOF group: %d ; in OTHER FOF groups: %d'%(sel.in_host_fof.sum(),(~sel.in_host_fof).sum()))
fof=np.array(mt.get_group_members(0,122)); 
print('  host FOF members total %d -> inside R200c %d, outside %d'%(len(fof),np.isin(sel.subhalo_id,fof).sum(),len(fof)-np.isin(sel.subhalo_id,fof).sum()))
lut=d.lookup; m=np.array([units.mass(d.SubhaloMass[lut[(122,int(i))]]) for i in sel.subhalo_id])
for lim in (1e7,1e8,1e9): print('  M_sub > %.0e Msun: %d  (of which other-FOF: %d)'%(lim,(m>lim).sum(),((m>lim)&~sel.in_host_fof).sum()))
o=~sel.in_host_fof
if o.any(): print('  other-FOF members: r/R200c range %.2f-%.2f; masses (Msun) max %.2e'%(sel.r_over_ref[o].min(),sel.r_over_ref[o].max(),m[o].max()))
# verify against the tracker's own inside_now flag for the FOF members
cfg=OrbitConfig(radius_key='R200c_from_M200'); tr=OrbitTracker(cosmo,units,cfg,boxsize=mt.BoxSize)
ids=sel.subhalo_id.tolist(); cat=tr.run([mt.get_track(int(i),122) for i in ids],host,ids=ids)
print('tracker inside_now over the selection: %d/%d inside (selection and tracker use identical z=0 geometry)'%(cat.scalar('inside_now').sum(),len(ids)))
print('  backsplash among them: %d ; infalling for first time: %d ; never crossed (born inside or tiny): n_infalls==0: %d'%(cat.scalar('is_backsplash').sum(),((cat.scalar('n_infalls')>=1)&cat.scalar('inside_now')).sum(),(cat.scalar('n_infalls')==0).sum()))
