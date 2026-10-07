#!/usr/bin/env python3
"""Diagnostic: history of the largest Mpeak/Mcurrent outlier and the frequency of single-snapshot peak spikes.

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/11_diag_mpeak_outlier.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import logging, numpy as np
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'05_plot_mpeak_vs_mcur.py')).read().split("print('inside R200c")[0])
i=int(np.argmax(np.where(cat.satellites('counts'),mpk/np.maximum(mcu,1),0)))
sid=ids[i]; t=mt.get_track(int(sid),122); print('outlier SubhaloNr',sid,' Mpeak %.2e Mcur %.2e'%(mpk[i],mcu[i]))
m=units.mass(t.Mass); L=t.extra['Len']
k=np.argmax(m); print('peak at snap',t.SnapNum[k],'a=%.3f'%t.Time[k],' Len at peak',L[k],' IsSubhalo at peak',bool(t.IsSubhalo[k]),' group M200 at peak %.2e'%units.mass(t.extra['GroupM200'][k]))
print('track length',len(t),'first snap',t.SnapNum[0])
print(' snap   a     Mbound     Len  IsSub  GroupNr')
for j in list(range(0,len(t),6))+[len(t)-1]:
    print(f' {t.SnapNum[j]:4d} {t.Time[j]:.3f} {m[j]:.2e} {L[j]:7d} {int(t.IsSubhalo[j])}  {t.extra["GroupNr"][j]}')
# how many others have a big single-snapshot mass spike?  ratio of peak to the median of neighbours
spk=[]
for j,s in enumerate(ids):
    if not cat.satellites('counts')[j]: continue
    tt=mt.get_track(int(s),122); mm=units.mass(tt.Mass); kk=np.argmax(mm)
    nb=np.r_[mm[max(kk-3,0):kk],mm[kk+1:kk+4]]
    spk.append(mm[kk]/np.median(nb) if nb.size else 1.0)
spk=np.array(spk); print('\npeak/median(+-3 snapshot neighbours): median %.2f; >2: %d; >5: %d of %d'%(np.median(spk),(spk>2).sum(),(spk>5).sum(),spk.size))
