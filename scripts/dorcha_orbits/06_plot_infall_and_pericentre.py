#!/usr/bin/env python3
"""Orbital parameters at first R200c crossing and pericentre vs current radius (R200c units, peculiar velocities; NOT in Benson's convention - see 08/09).

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/06_plot_infall_and_pericentre.py
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
from analysistools.orbits.cosmology import G_KPC_KMS2_MSUN as G
from analysistools.orbits.kinematics import host_reference
logging.basicConfig(level=logging.ERROR)
mt=MergerTreeTools(TREE,'SubFind',comoving=True,little_h=True,native_includes_h=True,native_is_comoving=True,loglevel=logging.ERROR)
cosmo=FlatLCDM.from_tree(mt)
units=TrackUnits(comoving=True,little_h=True,h=0.678,length_kpc=1000.0,mass_msun=1e10,velocity='physical')
host=mt.get_track(0,122)
sel=select_satellites_subfind(mt,host,cosmo,units,122); ids=sel.subhalo_id.tolist()
cfg=OrbitConfig(radius_key='R200c_from_M200',min_mass_fraction=1e-4,min_particles=50)
cat=OrbitTracker(cosmo,units,cfg,boxsize=mt.BoxSize).run([mt.get_track(int(i),122) for i in ids],host,ids=ids)
orb=cat.satellites('orbits'); cnt=cat.satellites('counts')
# host V200c(t), R200c(t)
th=cosmo.age(host.Time); R=host_reference(host,units,'R200c_from_M200',1.0,None,cosmo); M=units.mass(host.extra['GroupM200']); V=np.sqrt(G*M/R)
ti=cat.scalar('t_infall'); has=np.isfinite(ti)&np.isfinite(cat.scalar('vr_infall'))
vr=cat.scalar('vr_infall'); vt=cat.scalar('vt_infall'); L=np.linalg.norm(cat.scalar('L_infall'),axis=1)
V200=np.interp(np.nan_to_num(ti,nan=th[-1]),th,V); R200=np.interp(np.nan_to_num(ti,nan=th[-1]),th,R)
x_r=vr/V200; x_t=vt/V200; vtot=np.hypot(vr,vt)/V200; theta=np.degrees(np.arctan2(vt,np.abs(vr))); j=L/(R200*V200)
def stats(mask,name):
    k=mask&has; 
    print(f'{name}: N={k.sum()}  |v|/V200 median {np.median(vtot[k]):.2f} [{np.percentile(vtot[k],16):.2f},{np.percentile(vtot[k],84):.2f}] ; vr/V200 median {np.median(x_r[k]):+.2f} ; vt/V200 median {np.median(x_t[k]):.2f} ; theta median {np.median(theta[k]):.0f} deg ; j=L/(R200 V200) median {np.median(j[k]):.2f} [{np.percentile(j[k],16):.2f},{np.percentile(j[k],84):.2f}]; frac vr>0 (outgoing at crossing): {np.mean(x_r[k]>0):.2f}')
stats(orb,'orbit sample'); stats(cnt,'counts sample')
ko=orb&has; kc=cnt&has&~orb
fig,ax=plt.subplots(2,2,figsize=(11.5,8.6)); A=ax.ravel()
A[0].scatter(x_r[kc],x_t[kc],s=9,c='0.7',label=f'counts sample only (N={kc.sum()})',zorder=2)
A[0].scatter(x_r[ko],x_t[ko],s=26,c='crimson',edgecolor='k',lw=.3,label=f'orbit sample (N={ko.sum()})',zorder=3)
th_=np.linspace(0,np.pi,100)
for v in (0.5,1,1.5,2): A[0].plot(-v*np.cos(th_[::-1])*(-1)*-1 if False else v*np.cos(th_),v*np.sin(th_),'k:',lw=.6)
A[0].axvline(0,color='k',lw=.6); A[0].set_xlim(-2.6,1.6); A[0].set_ylim(0,2.2); A[0].set_aspect('equal')
A[0].set_xlabel('$v_r/V_{200c}$  (negative = infalling)'); A[0].set_ylabel('$v_t/V_{200c}$'); A[0].legend(fontsize=7.5,loc='upper left'); A[0].set_title('velocity at first $R_{200c}$ crossing (dotted: $|v|/V_{200c}$ = 0.5,1,1.5,2)',fontsize=9)
def hist(a,x,bins,xl,title=None):
    a.hist(x[kc],bins=bins,color='0.75',density=True,label='counts sample only')
    a.hist(x[ko],bins=bins,histtype='step',color='crimson',lw=1.8,density=True,label='orbit sample')
    a.axvline(np.median(x[ko]),color='crimson',ls='--',lw=1); a.axvline(np.median(x[kc|ko]),color='0.3',ls=':',lw=1)
    a.set_xlabel(xl); a.set_ylabel('probability density'); a.grid(alpha=.2,lw=.4)
hist(A[1],vtot,np.linspace(0,3,25),'$|v|/V_{200c}$'); A[1].legend(fontsize=7.5)
hist(A[2],theta,np.linspace(0,90,19),'infall angle $\\theta=\\arctan(v_t/|v_r|)$ [deg]  (0 = radial)')
hist(A[3],j,np.linspace(0,1.2,25),'$j=L/(R_{200c}V_{200c})$  (specific angular momentum at crossing)')
fig.suptitle('Orbital parameters at first $R_{200c}$ crossing: satellites inside $R_{200c}$ at $z=0$ (dashed red/dotted grey = medians)',fontsize=10)
fig.tight_layout(); fig.savefig(OUT+'infall_params.png',dpi=165)

# ---------- pericentre vs current radius
rn=cat.scalar('r_final'); rp_min=cat.scalar('r_peri_min'); rp_last=cat.scalar('r_peri_last'); nper=cat.scalar('n_pericentres')
R0=R[-1]
fig,ax=plt.subplots(1,2,figsize=(12.5,5.6),sharey=True)
for a_,rp,ttl in ((ax[0],rp_min,'deepest pericentre since infall'),(ax[1],rp_last,'most recent pericentre')):
    g=np.isfinite(rp)
    k1=g&cnt&~orb; k2=g&orb
    a_.scatter(rn[k1],rp[k1],s=10,c='0.7',zorder=2,label=f'counts sample only (N={k1.sum()})')
    sc=a_.scatter(rn[k2],rp[k2],s=34,c=np.minimum(nper[k2],4),cmap='viridis_r',vmin=1,vmax=4,edgecolor='k',lw=.35,zorder=3,label=f'orbit sample (N={k2.sum()})')
    x=np.geomspace(5,1000,10); a_.plot(x,x,'k-',lw=1,label='$r_{\\rm peri}=r_{\\rm now}$'); a_.plot(x,0.5*x,'k:',lw=.7); a_.plot(x,0.25*x,'k:',lw=.7)
    a_.axvline(R0,color='crimson',ls='--',lw=1,label='$R_{200c}$(z=0)')
    a_.set_xscale('log'); a_.set_yscale('log'); a_.set_xlim(8,600); a_.set_ylim(8,600); a_.set_aspect('equal'); a_.grid(alpha=.2,which='both',lw=.4)
    a_.set_xlabel('current radius, $r(z=0)$ [kpc]'); a_.set_title(ttl,fontsize=10)
    print(ttl,': orbit sample N=%d; r_now/r_peri median %.2f [16-84: %.2f,%.2f]; fraction currently within 1.2x of pericentre: %.2f; within 1.2x for counts sample: %.2f'%(k2.sum(),np.median(rn[k2]/rp[k2]),np.percentile(rn[k2]/rp[k2],16),np.percentile(rn[k2]/rp[k2],84),np.mean(rn[k2]/rp[k2]<1.2),np.mean(rn[g&cnt]/rp[g&cnt]<1.2)))
ax[0].set_ylabel('pericentric radius, $r_{\\rm peri}$ [kpc]'); ax[0].legend(fontsize=7.5,loc='upper left')
cb=fig.colorbar(sc,ax=ax,pad=0.01,shrink=.8); cb.set_label('number of pericentres (orbit sample; 4 = 4 or more)')
fig.suptitle('Pericentre vs current radius, satellites inside $R_{200c}$ at $z=0$ (dotted: $r_{\\rm peri}$ = 0.5, 0.25 $r_{\\rm now}$)',fontsize=10)
fig.savefig(OUT+'peri_vs_now.png',dpi=165)
nn=cnt&~np.isfinite(rp_min); print('counts-sample satellites with NO confirmed pericentre yet: %d of %d ; orbit sample: %d of %d'%(nn.sum(),cnt.sum(),(orb&~np.isfinite(rp_min)).sum(),orb.sum()))
print('orbit-sample satellites with r_now < r_peri_min (inconsistent?):',int(np.sum(orb&(rn<rp_min*0.98))))
