#!/usr/bin/env python3
"""Compare the four {r_vir,R200c} x {Hubble flow, none} variants with Benson's fit (needs 08 run for both modes).

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/09_benson_compare.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import pickle, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analysistools.orbits import benson as B
SC=OUT
E={m:pickle.load(open(SC+f'benson_events_{m}.pkl','rb')) for m in ('bn','200c')}
variants=[('bn',True,'$r_{vir}$ (BN98) + Hubble flow  [Benson definition]','k'),('bn',False,'$r_{vir}$, no Hubble flow','tab:blue'),('200c',True,'$R_{200c}$ + Hubble flow','tab:orange'),('200c',False,'$R_{200c}$, no Hubble flow  [requested]','crimson')]
wins={0.0:(0.0,0.25),0.5:(0.35,0.7),1.0:(0.8,1.3)}
def quant(ev,hub,m):
    Vr=ev['Vr'] if hub else ev['Vr_pec']; Vt=ev['Vt']; f2=ev['f2']
    phi=np.degrees(np.arctan2(Vt,Vr)); En=B.specific_energy(Vr,Vt,f2); e=B.eccentricity(Vr,np.maximum(Vt,1e-3),f2)
    return Vr,Vt,phi,En,e
def fit_stats(zf):
    v,_,pdf,_=B.benson05_grid(zf,3.0,601); w=pdf*(v[1]-v[0])**2; VR,VT=np.meshgrid(v,v,indexing='ij')
    En=B.specific_energy(VR,VT,1.0); en=B.eccentricity(VR,np.maximum(VT,1e-3),1.0); ph=np.degrees(np.arctan2(VT,VR))
    # median phi from the marginal
    _,_,_,pc,dph=B.benson05_marginals(zf); cdf=np.cumsum(dph)/dph.sum(); mphi=pc[np.searchsorted(cdf,0.5)]
    return (w*VR).sum(),(w*VT).sum(),mphi,w[En>0].sum(),w[np.abs(en-1)<0.1].sum()
print('Sample: first r_ref crossings, not a subhalo of a different FOF group, >=50 particles. Columns: N, <Vr>, <Vt>, median phi[deg], unbound frac, frac e in 1+-0.1')
for zf,(lo,hi) in wins.items():
    print('\n--- z in [%.2f,%.2f]  (compare Benson fit z=%.1f) ---'%(lo,hi,zf))
    for mode,hub,label,c in variants:
        ev=E[mode]; Vr,Vt,phi,En,e=quant(ev,hub,mode); m=ev['first']&ev['cen']&(ev['z']>=lo)&(ev['z']<hi)
        print('%-52s %4d  %.3f  %.3f  %5.1f   %.2f   %.2f'%(label.replace('$','').replace('\\\\',''),m.sum(),Vr[m].mean(),Vt[m].mean(),np.median(phi[m]),np.mean(En[m]>0),np.mean(np.abs(e[m]-1)<0.1)))
    r,t,mp,u,ep=fit_stats(zf); print('%-52s  --   %.3f  %.3f  %5.1f   %.2f   %.2f'%('Benson (2005) fit, f2=1',r,t,mp,u,ep))
fig,ax=plt.subplots(2,3,figsize=(13.5,7.6))
for row,zf in enumerate((0.0,1.0)):
    lo,hi=wins[zf]; v,dvr,dvt,pc,dph=B.benson05_marginals(zf)
    for col,(key,bins,fx,fy,xl) in enumerate((('Vr',np.linspace(-0.2,2.2,25),v,dvr,'$V_r$ (inward) / $V_{ref}$'),('Vt',np.linspace(0,2.2,23),v,dvt,'$V_\\theta$ / $V_{ref}$'),('phi',np.linspace(0,90,19),pc,dph,'infall angle $\\phi$ [deg]'))):
        a_=ax[row,col]
        for mode,hub,label,c in variants:
            ev=E[mode]; Vr,Vt,phi,En,e=quant(ev,hub,mode); m=ev['first']&ev['cen']&(ev['z']>=lo)&(ev['z']<hi)
            x={'Vr':Vr,'Vt':Vt,'phi':phi}[key]
            a_.hist(x[m],bins=bins,density=True,histtype='step',color=c,lw=2.2 if mode=='200c' and not hub else 1.3,label=label if (row==0 and col==0) else None)
        a_.plot(fx,fy,'-',color='0.35',lw=3,alpha=.6,label='Benson (2005) fit' if (row==0 and col==0) else None)
        a_.set_xlabel(xl); a_.grid(alpha=.2,lw=.4)
        if col==0: a_.set_ylabel('probability density')
    ax[row,0].set_title(f'z in [{lo:g},{hi:g}] vs fit z={zf:g}',fontsize=9,loc='left')
ax[0,0].legend(fontsize=7,loc='upper left')
fig.suptitle('Dorcha host 0: effect of Hubble flow and radius definition on infall parameters (first crossings, not pre-processed)',fontsize=10)
fig.tight_layout(); fig.savefig(SC+'benson_variants.png',dpi=160)
