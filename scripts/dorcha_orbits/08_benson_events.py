#!/usr/bin/env python3
"""All first crossings of r_ref by the host's subhaloes over cosmic history, in Benson (2005) quantities. Usage: 08_benson_events.py bn|200c  (bn = Bryan-Norman r_vir, 200c = R200c). Writes benson_events_<mode>.pkl.

Run from the repository root:  PYTHONPATH=. python scripts/dorcha_orbits/08_benson_events.py
See scripts/dorcha_orbits/README.md for the environment variables.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import TREE, OUT, OUT_DIR
import logging, time, sys, numpy as np, h5py, pickle
MODE=sys.argv[1] if len(sys.argv)>1 else 'bn'   # 'bn' = Bryan-Norman r_vir ; '200c' = R200c
from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.orbits import *
from analysistools.orbits import benson as B
from analysistools.orbits.cosmology import G_KPC_KMS2_MSUN as G
logging.basicConfig(level=logging.ERROR)
SC=OUT
TREE=TREE
mt=MergerTreeTools(TREE,'SubFind',comoving=True,little_h=True,native_includes_h=True,native_is_comoving=True,loglevel=logging.ERROR)
cosmo=FlatLCDM.from_tree(mt); d=mt.data; h=0.678
units=TrackUnits(comoving=True,little_h=True,h=h,length_kpc=1000.0,mass_msun=1e10,velocity='physical')
host=mt.get_track(0,122)
a=np.asarray(host.Time); z=1/a-1; hs=np.asarray(host.SnapNum); k_of={int(s):i for i,s in enumerate(hs)}
# ---- host virial radius in Benson's sense (BN98 top-hat), from NFW(M200c, R200c, rs=Rmax/2.163)
R200=host_reference(host,units,'R200c_from_M200',1.0,None,cosmo) if False else None
from analysistools.orbits.kinematics import host_reference
R200=host_reference(host,units,'R200c_from_M200',1.0,None,cosmo)
M200=units.mass(host.extra['GroupM200']); rmax=units.physical_kpc(host.extra['VmaxRad'],a)
ok=(M200>0)&(rmax>0)
delta=B.bryan_norman_delta(cosmo,z) if MODE=='bn' else np.full(a.size,200.0)
rvir=np.full(a.size,np.nan)
if MODE=='bn': rvir[ok]=B.nfw_virial_radius(M200[ok],R200[ok],rmax[ok]/2.1626,delta[ok],cosmo,a[ok])
else: rvir[ok]=R200[ok]
print('z=0: Delta_vir=%.1f  R200c=%.1f kpc  r_vir=%.1f kpc  (ratio %.3f)'%(delta[-1],R200[-1],rvir[-1],rvir[-1]/R200[-1]))
for c in ((5,10,20) if MODE=='bn' else ()):                       # sensitivity of r_vir/R200c to the unmeasured concentration
    rv=B.nfw_virial_radius(M200[-1],R200[-1],R200[-1]/c,delta[-1],cosmo,1.0)[0]; print('   if c200=%2d -> r_vir/R200c=%.3f'%(c,rv/R200[-1]))
host.extra['Rvir_BN']=np.where(np.isfinite(rvir),rvir,0.0)*h/(1000.0*a)       # tree units (comoving Mpc/h)
# ---- candidate first-crossing events at snapshot level
NMIN=50
snap=np.asarray(d.SnapNum); pos=np.asarray(d.SubhaloPos,float); prog=np.asarray(d.TreeMainProgenitor); L=np.asarray(d.SubhaloLen); sn=np.asarray(d.SubhaloNr)
order=np.argsort(snap,kind='stable'); bounds=np.searchsorted(snap[order],np.arange(0,124))
box=mt.BoxSize
def rphys(p,s):                                   # physical kpc from host centre at snapshot s
    k=k_of[s]; dx=p-host.Pos[k]; dx=(dx+0.5*box)%box-0.5*box
    return np.linalg.norm(dx,axis=1)*a[k]*1000.0/h
cands=[]
for s in range(hs.min()+1,123):
    k=k_of[s]
    if z[k]>1.45 or not np.isfinite(rvir[k]) or not np.isfinite(rvir[k_of[s-1]]): continue
    rows=order[bounds[s]:bounds[s+1]]; pr=prog[rows]; has=pr>=0
    rows,pr=rows[has],pr[has]; good=snap[pr]==s-1; rows,pr=rows[good],pr[good]
    rs_=rphys(pos[rows],s); rp_=rphys(pos[pr],s-1)
    m=(rs_<rvir[k])&(rp_>=rvir[k_of[s-1]])&(L[pr]>=NMIN)&(sn[rows]!=host.extra['SubhaloNr'][k])
    cands+= [(int(sn[r]),s,int(r)) for r in rows[m]]
print('candidate first-crossing nodes (>=%d particles before crossing, z<1.45): %d'%(NMIN,len(cands)))
# ---- survivors: does the branch reach z=0, and is it inside r_vir there?
desc=h5py.File(TREE,'r')['TreeHalos/TreeDescendant'][()]
rows=np.array([c[2] for c in cands]); cur=rows.copy()
for _ in range(125):
    nxt=np.where(desc[cur]>=0,desc[cur],cur); 
    if np.array_equal(nxt,cur): break
    cur=nxt
reach0=snap[cur]==122
r0=np.full(len(cands),np.inf); r0[reach0]=rphys(pos[cur[reach0]],122)
surv=reach0&(r0<rvir[-1])
# ---- orbits at crossing
t0=time.time()
tracks=[mt.get_track(c[0],c[1]) for c in cands]
cfg=OrbitConfig(radius_key='Rvir_BN',hysteresis=0.01,orphans=False)
cat=OrbitTracker(cosmo,units,cfg,boxsize=mt.BoxSize).run(tracks,host,ids=list(range(len(cands))))
print('tracked in %.0fs'%(time.time()-t0))
tin=cat.scalar('t_infall'); vr=cat.scalar('vr_infall'); vt=cat.scalar('vt_infall'); rv=cat.scalar('r_ref_infall'); msat=cat.scalar('mass_infall')
sc_snap=np.array([c[1] for c in cands]); t_s=cosmo.age(a[[k_of[s] for s in sc_snap]]); t_p=cosmo.age(a[[k_of[s-1] for s in sc_snap]])
first=np.isfinite(tin)&(tin>t_p-1e-9)&(tin<=t_s+1e-9)&np.isfinite(vr)&~cat.scalar('born_inside')
print('events whose FIRST crossing falls in the candidate interval:',first.sum(),' (dropped: re-entries / no crossing:',(~first).sum(),')')
ac=cosmo.a_of_t(np.where(first,tin,1.0)); zc=1/ac-1
Dc=B.bryan_norman_delta(cosmo,zc) if MODE=='bn' else np.full_like(zc,200.0); Vvir=B.v_vir_from_r_vir(cosmo,ac,rv,Dc)
Hc=cosmo.hubble(ac)/1.0227121650537077
Vr=-(vr+Hc*rv)/Vvir; Vt=vt/Vvir; Vr_pec=-vr/Vvir
Mvir=Vvir**2*rv/G; f2=B.f2_of(msat,Mvir)
cen_strict=np.array([not t.IsSubhalo[-2] if len(t)>1 else True for t in tracks])
# not a subhalo of a DIFFERENT FOF group (being in the host's own group before crossing is fine)
other=np.array([bool(t.IsSubhalo[-2] and t.HostID[-2]!=host.extra['SubhaloNr'][k_of[int(t.SnapNum[-2])]]) if len(t)>1 else False for t in tracks])
cen=~other
ev=dict(cen_strict=cen_strict,first=first,z=zc,Vr=Vr,Vt=Vt,Vr_pec=Vr_pec,f2=f2,cen=cen,surv=surv,snap=sc_snap,id=np.array([c[0] for c in cands]),msat=msat,Vvir=Vvir,Hc=Hc,rv=rv)
pickle.dump(ev,open(SC+'benson_events_%s.pkl'%MODE,'wb'))
g=first
print('events: N=%d ; central at infall: %d ; survivors-at-z0: %d'%(g.sum(),(g&cen).sum(),(g&surv).sum()))
print('Hubble-flow shift on V_r: mean %.3f V_vir (=sqrt(2/Delta)) ; with vs without: median Vr %.2f vs %.2f'%(np.mean(Hc[g]*rv[g]/Vvir[g]),np.median(Vr[g]),np.median(Vr_pec[g])))
