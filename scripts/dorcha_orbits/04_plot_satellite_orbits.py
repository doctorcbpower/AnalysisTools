#!/usr/bin/env python3
"""Radius vs expansion factor of the ten most massive satellites (inside R200c at z=0) with host R200c(a). Usage: 04_plot_satellite_orbits.py [a_min a_max] (default 0.5 1.0).

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/04_plot_satellite_orbits.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import logging, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.orbits import *
from analysistools.orbits.kinematics import host_reference
logging.basicConfig(level=logging.ERROR)
A_MIN=float(sys.argv[1]) if len(sys.argv)>2 else 0.5
A_MAX=float(sys.argv[2]) if len(sys.argv)>2 else 1.0
YMAX=[0.0]
mt=MergerTreeTools(TREE,'SubFind',comoving=True,little_h=True,native_includes_h=True,native_is_comoving=True,loglevel=logging.ERROR)
cosmo=FlatLCDM.from_tree(mt); d=mt.data; lut=d.lookup
units=TrackUnits(comoving=True,little_h=True,h=0.678,length_kpc=1000.0,mass_msun=1e10,velocity='physical')
host=mt.get_track(0,122)
sel=select_satellites_subfind(mt,host,cosmo,units,122)               # inside R200c, any FOF group
m=np.array([units.mass(d.SubhaloMass[lut[(122,int(i))]]) for i in sel.subhalo_id])
order=np.argsort(m)[::-1]; top=sel.subhalo_id[order[:10]]
tracks=[mt.get_track(int(i),122) for i in top]
cfg=OrbitConfig(radius_key='R200c_from_M200',min_mass_fraction=1e-4,min_particles=50)
cat=OrbitTracker(cosmo,units,cfg,boxsize=mt.BoxSize).run(tracks,host,ids=top.tolist())
rref=host_reference(host,units,'R200c_from_M200',1.0,None,cosmo)
cols=plt.get_cmap('tab10').colors
fig,ax=plt.subplots(figsize=(11,5.8))
for i,(sid,tr) in enumerate(zip(top,tracks)):
    ser=host_frame_series(tr,host,cosmo,units,mt.BoxSize,'R200c_from_M200'); fine=interpolate_orbit(ser,cosmo,16)
    ax.plot(fine.a,fine.r,'-',color=cols[i],lw=1.3)
    w_=(fine.a>=A_MIN)&(fine.a<=A_MAX)
    if w_.any(): YMAX[0]=max(YMAX[0],float(fine.r[w_].max()))
    ax.plot(ser.a,ser.r,'o',color=cols[i],ms=3,mec='none',label=f"{int(sid)}  ({m[order[i]]:.1e} M$_\\odot$)")
    t_inf=cat.record(i).scalars['t_infall']
    for e in cat.record(i).events:
        if e.a<A_MIN or e.a>A_MAX: continue
        if e.kind==PERICENTRE and not (np.isfinite(t_inf) and e.t>=t_inf): continue   # only orbits about the host
        mk={PERICENTRE:'D',CROSS_IN:'v',CROSS_OUT:'^'}.get(e.kind)
        if mk: ax.plot(e.a,e.r,mk,color=cols[i],ms=7 if mk=='D' else 8,mec='k',mew=.8,zorder=6)
ax.plot(host.Time,rref,'k--',lw=2.6,label='host $R_{200c}$',zorder=5)
ax.set_xlim(A_MIN,A_MAX); ax.set_ylim(0,1.05*YMAX[0])
ax.set_xlabel('expansion factor, $a$'); ax.set_ylabel('halo-centric radius, $r$ [physical kpc]'); ax.grid(alpha=.25,lw=.5)
ax2=ax.secondary_xaxis('top',functions=(lambda a:1/a-1,lambda z:1/(1+z))); zt=[zz for zz in (10,5,3,2,1.5,1,0.7,0.5,0.4,0.3,0.2,0.1,0.05,0) if A_MIN<=1/(1+zz)<=A_MAX]; ax2.set_xticks(zt); ax2.set_xticklabels([f'{zz:g}' for zz in zt]); ax2.set_xlabel('redshift, $z$')
ax.set_title(f'Ten most massive satellites (inside $R_{{200c}}$ at $z=0$) of Dorcha host 0, $a$={A_MIN:g}-{A_MAX:g}',fontsize=11)
leg=ax.legend(fontsize=7.5,loc='upper left',bbox_to_anchor=(1.02,1.0),title='SubhaloNr (mass at z=0)',title_fontsize=7.5,framealpha=.92); ax.add_artist(leg)
ax.legend(handles=[Line2D([],[],marker='o',color='gray',ls='',ms=4,label='tree snapshot'),Line2D([],[],color='gray',label='Hermite-interpolated orbit'),
 Line2D([],[],marker='D',color='w',mfc='gray',mec='k',ls='',label='pericentre (after first infall)'),Line2D([],[],marker='v',color='w',mfc='gray',mec='k',ls='',ms=8,label='crosses in through $R_{200c}$'),
 Line2D([],[],marker='^',color='w',mfc='gray',mec='k',ls='',ms=8,label='crosses out')],fontsize=7.5,loc='lower left',bbox_to_anchor=(1.02,0.0),framealpha=.92)
fig.tight_layout(); out=os.path.join(OUT_DIR,f'satellite_orbits_a{A_MIN:g}_{A_MAX:g}.png'); fig.savefig(out,dpi=170); print(out)
print('top10 (inside R200c):',[int(i) for i in top])
for i,sid in enumerate(top):
    sc=cat.record(i).scalars
    print(int(sid),'M=%.1e'%m[order[i]],'r_z0=%.0f (%.2f R200c)'%(sc['r_final'],sc['r_final']/rref[-1]),'infall a=%s'%(round(float(cosmo.a_of_t(sc['t_infall'])),3) if np.isfinite(sc['t_infall']) else None),'n_peri',int(sc['n_pericentres']),'r_peri_min',None if not np.isfinite(sc['r_peri_min']) else round(float(sc['r_peri_min']),0),'preproc',bool(sc['preprocessed']),'backspl-before' ,int(sc['n_outfalls']))
# new default samples on the real tree (everyone inside R200c)
ids=sel.subhalo_id.tolist(); full=OrbitTracker(cosmo,units,cfg,boxsize=mt.BoxSize).run([mt.get_track(int(i),122) for i in ids],host,ids=ids)
print('\nhost M200c(z=0)=%.3e Msun ; effective mass cut = %.3e Msun (= %.0f particles; 50-particle cut = %.2e Msun)'%(full.metadata['host_m200_final_msun'],full.metadata['min_mass_effective_msun'],full.metadata['min_mass_effective_msun']/2.324e5,50*2.324e5))
print('satellites inside R200c: %d ; counts sample (peak basis): %d ; orbit sample (final basis): %d'%(len(ids),full.satellites('counts').sum(),full.satellites('orbits').sum()))
print('only the particle cut (50): counts %d, orbits %d'%(full.select(min_particles=50,basis='peak',inside_only=True).sum(),full.select(min_particles=50,basis='final',inside_only=True).sum()))
