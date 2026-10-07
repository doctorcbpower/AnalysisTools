#!/usr/bin/env python3
"""Satellite counts for user-defined particle/mass thresholds, z=0 basis vs peak/final catalogue basis.

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/07b_threshold_counts.py
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
mp=np.median(np.asarray(d.SubhaloMass)[np.asarray(d.SubhaloLen)>200]/np.asarray(d.SubhaloLen)[np.asarray(d.SubhaloLen)>200])
print('bound mass per particle (tree units): %.3e  -> %.3e Msun (incl. 1/h conversion)'%(mp,units.mass(mp)))
base=select_satellites_subfind(mt,host,cosmo,units,122); print('no threshold:',len(base))
print('\nz=0 basis, applied at selection time:')
for npn in (None,20,32,50,100,300):
    print('  min_particles=%-5s -> %4d'%(npn,len(select_satellites_subfind(mt,host,cosmo,units,122,min_particles=npn))))
for mm in (1e7,1e8,1e9):
    print('  min_mass=%.0e Msun -> %4d'%(mm,len(select_satellites_subfind(mt,host,cosmo,units,122,min_mass=mm))))
ids=base.subhalo_id.tolist()
cfg=OrbitConfig(radius_key='R200c_from_M200'); cat=OrbitTracker(cosmo,units,cfg,boxsize=mt.BoxSize).run([mt.get_track(int(i),122) for i in ids],host,ids=ids)
print('\ncatalogue mask, same satellites (inside R200c), PEAK vs FINAL basis:')
print('  %-18s %8s %8s'%('threshold','peak','final'))
for npn in (20,32,50,100,300):
    print('  min_particles=%-4d %8d %8d'%(npn,cat.select(min_particles=npn,basis='peak').sum(),cat.select(min_particles=npn,basis='final').sum()))
for mm in (1e7,1e8,1e9):
    print('  min_mass=%.0e     %8d %8d'%(mm,cat.select(min_mass=mm,basis='peak').sum(),cat.select(min_mass=mm,basis='final').sum()))
npk=cat.scalar('npart_peak'); npf=cat.scalar('npart_final'); print('\nmedian peak/final particle counts:',np.median(npk),np.median(npf),' ; fraction with peak>=2x final:',np.mean(npk>=2*npf).round(2))
