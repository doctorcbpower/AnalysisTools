#!/usr/bin/env python3
"""Mpeak vs Mcurrent for the counts sample (Mpeak >= resolution cut), plus Mcur/Mpeak histogram.

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/05_plot_mpeak_vs_mcur.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import logging, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.orbits import *
logging.basicConfig(level=logging.ERROR)
mt=MergerTreeTools(TREE,'SubFind',comoving=True,little_h=True,native_includes_h=True,native_is_comoving=True,loglevel=logging.ERROR)
cosmo=FlatLCDM.from_tree(mt)
units=TrackUnits(comoving=True,little_h=True,h=0.678,length_kpc=1000.0,mass_msun=1e10,velocity='physical')
host=mt.get_track(0,122)
sel=select_satellites_subfind(mt,host,cosmo,units,122)
ids=sel.subhalo_id.tolist()
cfg=OrbitConfig(radius_key='R200c_from_M200',min_mass_fraction=1e-4,min_particles=50)
cat=OrbitTracker(cosmo,units,cfg,boxsize=mt.BoxSize).run([mt.get_track(int(i),122) for i in ids],host,ids=ids)
cut=cat.metadata['min_mass_effective_msun']
mpk=cat.scalar('mass_peak'); mcu=cat.scalar('mass_final'); npk=cat.scalar('npart_peak')
inside=cat.scalar('inside_now'); cnt=cat.satellites('counts'); orb=cat.satellites('orbits')
a_in=np.where(np.isfinite(cat.scalar('t_infall')),cosmo.a_of_t(np.nan_to_num(cat.scalar('t_infall'),nan=1.0)),np.nan)
print('inside R200c:',inside.sum(),' counts sample (Mpeak>=%.2e):'%cut,cnt.sum(),' also Mcur>=cut (orbit sample):',orb.sum())
assert (mpk[cnt]>=cut).all() and (mcu<=mpk*(1+1e-9)).all()
f=mcu[cnt]/mpk[cnt]
print('Mcur/Mpeak over counts sample: median %.3f, 16-84%% [%.3f, %.3f], min %.4f; fraction with Mcur/Mpeak<0.1: %.2f; <0.01: %.2f'%(np.median(f),*np.percentile(f,[16,84]),f.min(),np.mean(f<0.1),np.mean(f<0.01)))
print('fraction of counts sample dropping below the cut by z=0: %.2f'%np.mean(~orb[cnt]))
print('Mpeak range %.2e-%.2e ; Mcur range %.2e-%.2e ; Mcur<32 particles equiv(7.4e6): %d'%(mpk[cnt].min(),mpk[cnt].max(),mcu[cnt].min(),mcu[cnt].max(),(mcu[cnt]<32*2.324e5).sum()))
# Mcur/Mpeak vs infall time
for lo,hi in [(0.0,0.5),(0.5,0.7),(0.7,0.85),(0.85,1.01)]:
    k=cnt&(a_in>=lo)&(a_in<hi); 
    if k.sum(): print('  infall a in [%.2f,%.2f): N=%3d  median Mcur/Mpeak=%.3f'%(lo,hi,k.sum(),np.median(mcu[k]/mpk[k])))
k=cnt&~np.isfinite(a_in); print('  no recorded infall crossing: N=%d, median Mcur/Mpeak=%.3f'%(k.sum(),np.median(mcu[k]/mpk[k]) if k.sum() else np.nan))

fig,ax=plt.subplots(1,2,figsize=(12.5,5.4),gridspec_kw={'width_ratios':[1.15,1]})
a0=ax[0]; lo=1e7
kk=cnt&np.isfinite(a_in); sc=a0.scatter(mpk[kk],mcu[kk],c=a_in[kk],cmap='viridis_r',s=22,vmin=0.3,vmax=1.0,edgecolor='k',linewidth=.25,zorder=3)
k2=cnt&~np.isfinite(a_in); a0.scatter(mpk[k2],mcu[k2],c='0.6',s=22,edgecolor='k',linewidth=.25,zorder=3,label='no recorded R$_{200c}$ crossing')
ring=cnt&orb; a0.scatter(mpk[ring],mcu[ring],s=70,facecolor='none',edgecolor='crimson',linewidth=.9,zorder=4,label=f'also $M_{{\\rm cur}}\\geq$ cut ({orb.sum()}: orbit sample)')
x=np.geomspace(1e8,1.5e12,10); a0.plot(x,x,'k-',lw=1,label='$M_{\\rm cur}=M_{\\rm peak}$'); a0.plot(x,0.1*x,'k:',lw=.8); a0.plot(x,0.01*x,'k:',lw=.8)
a0.text(2.5e11,2.5e10,'0.1',fontsize=8,rotation=33); a0.text(2.5e11,2.5e9,'0.01',fontsize=8,rotation=33)
j=int(np.argmax(np.where(cnt,mpk/np.maximum(mcu,1),0)))
a0.annotate("SubhaloNr %d: ~1:2 merger remnant\n(central of its own FOF group until a=0.67)"%ids[j],xy=(mpk[j],mcu[j]),xytext=(1.5e10,1.6e6),fontsize=7.5,arrowprops=dict(arrowstyle="->",lw=.8))
a0.axvline(cut,color='crimson',ls='--',lw=1); a0.axhline(cut,color='crimson',ls='--',lw=1)
a0.set_xscale('log'); a0.set_yscale('log'); a0.set_xlim(1.5e8,1.5e12); a0.set_ylim(1e6,5e11)
a0.set_xlabel('$M_{\\rm peak}$ [M$_\\odot$] (peak bound subhalo mass along the main branch)'); a0.set_ylabel('$M_{\\rm current}$ [M$_\\odot$] (bound mass at $z=0$)')
cb=fig.colorbar(sc,ax=a0,pad=0.01); cb.set_label('expansion factor at first $R_{200c}$ crossing')
a0.legend(fontsize=7.5,loc='upper left',framealpha=.92); a0.grid(alpha=.2,which='both',lw=.4)
a0.set_title(f'Satellites inside $R_{{200c}}$ with $M_{{\\rm peak}}\\geq${cut:.2g} M$_\\odot$ (N={cnt.sum()})',fontsize=10)
a1=ax[1]; bins=np.geomspace(max(f.min(),1e-4),1.0,22)
a1.hist(f,bins=bins,color='0.55',edgecolor='k',lw=.4); a1.set_xscale('log'); a1.axvline(np.median(f),color='crimson',lw=1.5,label=f'median = {np.median(f):.2f}')
a1.set_xlabel('$M_{\\rm current}/M_{\\rm peak}$'); a1.set_ylabel('number of satellites'); a1.legend(fontsize=8); a1.grid(alpha=.2,lw=.4)
a1.set_title('Surviving bound-mass fraction',fontsize=10)
fig.tight_layout(); out=os.path.join(OUT_DIR,'mpeak_vs_mcur.png'); fig.savefig(out,dpi=170); print(out)
