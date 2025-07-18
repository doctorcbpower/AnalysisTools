#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:45:15 2023

@author: 00075868
"""

import h5py
import numpy as np
import os
import snapshot_tools as st

class ProfileTools:
    def __init__(self,**kwargs):
        self.comoving_units=True
        # Use comoving or physical units when plotting?
        if 'comoving_units' in kwargs:    
            self.comoving_units=kwargs.get('comoving_units')
        self.fr_cut=5             
        # Size of region in units of R200
        if 'fr_cut' in kwargs:
            self.radius_cut=kwargs.get('fr_cut')
        self.numbins=25
        # Number of bins to plot
        if 'numbins' in kwargs:
            self.numbins=kwargs.get('numbins')
        self.halo_id=0
        # Which halo to plot?
        if 'halo_id' in kwargs:
            self.halo_id=kwargs.get('halo_id')
            
    ''' Load particles in a region defined by the properties of 
        the halo in the halo catalogue. 
        pos_centre=3D vector defining centre
        size=3D vector defining lengths along x,y,z directions
        region_type=sphere/cube/FOF - spherical region, cuboidal
        region (with dimensions set by user), FOF only
    '''
        
    def LoadParticlesInRegion(self,particle_data,pos_centre,size,region_type='sphere',**kwargs):
        if region_type=='FOF':
            if 'halo_data' in kwargs:
                halo_data=kwargs.get('halo_data')
                ikeep=np.isin(particle_data.pids,halo_data.particle_id)
            else:
                print('Need to pass halo_data')
        else:
            delta=particle_data.pos-pos_centre
            if region_type=='cuboid':
                ikeep=np.logical_and(np.abs(delta[:,0])<=self.fr_cut*size[0],np.abs(delta[:,1])<=self.fr_cut*size[1])
                ikeep=np.logical_and(ikeep,np.abs(delta[:,2])<=self.fr_cut*size[2])            
            if region_type=='sphere':
                r2=delta[:,0]**2+delta[:,1]**2+delta[:,2]**2
                ikeep=np.where(r2<=(self.fr_cut*size[0])**2)[0]
        self.keepids=particle_data.pids[ikeep]

    def SphericallyBinParticlesInRegion(self,particle_data,ikeep,pos_centre,rlim,**kwargs):
        self.r=np.sqrt(particle_data.pos[ikeep][:,0]**2+particle_data.pos[ikeep][:,1]**2+particle_data.pos[ikeep][:,2]**2)
        self.lrmin=np.log10(rlim[0])
        self.lrmax=np.log10(rlim[1])
        ibins=np.linspace(self.lrmin,self.lrmax,self.numbins)
        self.ilist=np.digitize(np.log10(self.r),ibins)
        
#    def select_particles(self,val,valoffset,size,geometry,**kwargs):
#        dval=val-valoffset
#        # First check for periodicity
#        if kwargs.get('periodic')==True and kwargs.get('scale_length')!=None:
#            scale_length=kwargs.get('scale_length')
#            dval=np.where(dval>0.5*scale_length,dval-scale_length,dval)
#            dval=np.where(dval<-0.5*scale_length,dval+scale_length,dval)
#        else:
#            print('Ignoring periodicity')
#        # Impose cut based on desired geometry
#        if geometry=='cubic':
#            ipick=np.logical_and(np.abs(dval[:,0])<size,np.abs(dval[:,1])<size)
#            ipick=np.logical_and(ipick,np.abs(dval[:,2])<size)
#        elif geometry=='spherical':
#            r2=dval[:,0]**2+dval[:,1]**2+dval[:,2]**2
#            ipick=np.where(r2<size*size)[0]
#        else:
#            print('Undefined geometry')
#            ipick=None
#        return ipick
#
    def radial_bins(self,r,rlim,nbins,bintype='log'):
        if bintype=='log':
            lrmin=np.log10(rlim[0])
            lrmax=np.log10(rlim[1])
            lrbins=np.linspace(lrmin,lrmax,nbins)
            rbins=10**lrbins
        else:
            rbins=np.linspace(rlim[0],rlim[1],nbins)
        indx=np.digitize(r,bins=rbins,right=False)

        return indx,rbins

    def radial_profile(self,r,nbins,rbins,index,var,**kwargs):
        r_av=np.zeros(nbins,dtype=np.float64)    # Average value of radius in bin
        var_av=np.zeros(nbins,dtype=np.float64)  # Average value of variable in bin
        var_med=np.zeros(nbins,dtype=np.float64) # Median value of variable in bin
        var_percentiles=np.zeros([nbins,2],dtype=np.float64)  # Upper and lower percentiles in bin

        for i in range(nbins):
            r_av[i]=np.mean(r[index==i])
            var_sum=np.sum(var[index==i])
            if i==0:
                volume=(4*np.pi/3.)*rbins[i]**3
            else:
                volume=(4*np.pi/3.)*(rbins[i]**3-rbins[i-1]**3)
                var_av[i]=var_sum/volume

            var_med[i]=np.median(var[index==i])
        
        return r_av,var_av,var_med

    def plot_mass_density_profile(self,pos,mass,poscen,size,lbox,rlim,numbins,geometry='spherical',type='average'):
        ipick=st.select_particles(pos,poscen,size,geometry,periodic=True,scale_length=lbox)
        dpos=pos[ipick]-poscen
        dpos=np.where(dpos>0.5*lbox,dpos-lbox,dpos)
        dpos=np.where(dpos<-0.5*lbox,dpos+lbox,dpos)
        r=np.sqrt(dpos[:,0]**2+dpos[:,1]**2+dpos[:,2]**2)
<<<<<<< HEAD
        indx,rbins=self.radial_bins(r=r,rmin=rlim[0],rmax=rlim[1],nbins=numbins,bintype='log')
=======
        indx,rbins=self.radial_bins(r=r,rlim=rlim,nbins=numbins,bintype='log')
>>>>>>> 5bb954d (Synchronising after too long - will start to tidy up scripts and document them)
        rav,rhoav,rhomed=self.radial_profile(r=r,nbins=numbins,rbins=rbins,index=indx,var=mass[ipick])
        if type=='average':
            return rav,rhoav
        else:
            return rav,rhomed
        
<<<<<<< HEAD
    def kinematic_radial_profile(self,pos,vel,poscen,velcen,size,lbox,rlim,numbins,bintype='log',**kwargs):
=======
    def plot_gas_profile(self,pos,rho,u,ipick,poscen,size,lbox,rlim,numbins,geometry='spherical',type='average'):
#        ipick=st.select_particles(pos,poscen,size,geometry,periodic=True,scale_length=lbox)
        dpos=pos[ipick]-poscen
        dpos=np.where(dpos>0.5*lbox,dpos-lbox,dpos)
        dpos=np.where(dpos<-0.5*lbox,dpos+lbox,dpos)
        r=np.sqrt(dpos[:,0]**2+dpos[:,1]**2+dpos[:,2]**2)
        indx,rbins=self.radial_bins(r=r,rlim=rlim,nbins=numbins,bintype='linear')
        rav,rhoav,rhomed=self.radial_profile(r=r,nbins=numbins,rbins=rbins,index=indx,var=rho[ipick])
        ruav,uav,umed=self.radial_profile(r=r,nbins=numbins,rbins=rbins,index=indx,var=u[ipick])
        if type=='average':
            return rav,rhoav,ruav,uav
        else:
            return rav,rhomed,ruav,umed
        
    def plot_kinematic_profile(self,pos,vel,poscen,velcen,size,lbox,rlim,numbins,bintype='log'):
        print(pos,vel,poscen,velcen,size,lbox,rlim,numbins,bintype)
>>>>>>> 5bb954d (Synchronising after too long - will start to tidy up scripts and document them)
        ipick=st.select_particles(pos,poscen,size,geometry='spherical',periodic=True,scale_length=lbox)
        dpos=pos[ipick]-poscen
        dpos=np.where(dpos>0.5*lbox,dpos-lbox,dpos)
        dpos=np.where(dpos<-0.5*lbox,dpos+lbox,dpos)
        r=np.sqrt(dpos[:,0]**2+dpos[:,1]**2+dpos[:,2]**2)
        dvel=vel[ipick]-velcen
        
        if bintype=='log':
            lrmin=np.log10(rlim[0])
            lrmax=np.log10(rlim[1])
            lrbins=np.linspace(lrmin,lrmax,numbins)
            rbins=10**lrbins
        else:
            rbins=np.linspace(rlim[0],rlim[1],numbins)
        indx=np.digitize(r,bins=rbins,right=False)

        rav=np.zeros(numbins,dtype=np.float64)
        vrav=np.zeros(numbins,dtype=np.float64)
        sigmar_av=np.zeros(numbins,dtype=np.float64)
        sigma_av=np.zeros(numbins,dtype=np.float64)
    
        for i in range(numbins):
            rav[i]=np.mean(r[indx==i])
            vr=dpos[:,0][indx==i]*dvel[:,0][indx==i]+dpos[:,1][indx==i]*dvel[:,1][indx==i]+dpos[:,2][indx==i]*dvel[:,2][indx==i]
            vr=vr/r[indx==i]
            vr2=vr*vr
            vrav[i]=np.mean(vr)
            vr2av=np.mean(vr*vr)
            sigmar_av[i]=np.sqrt(vr2av-vrav[i]**2)
            vxav=np.mean(dvel[:,0][indx==i])
            vyav=np.mean(dvel[:,1][indx==i])
            vzav=np.mean(dvel[:,2][indx==i])
            vx2av=np.mean(dvel[:,0][indx==i]**2)
            vy2av=np.mean(dvel[:,1][indx==i]**2)
            vz2av=np.mean(dvel[:,2][indx==i]**2)
            sigma_av[i]=np.sqrt((vx2av-vxav*vxav)+(vy2av-vyav*vyav)+(vz2av-vzav*vzav))
        return rav,vrav,sigmar_av,sigma_av

# rmed_mr=np.array([])
# umed_mr=np.array([])
# rhomed_mr=np.array([])
# rhoshell_mr=np.array([])

# for i in range(numbins):
#     if np.any(ilist==i)==True:
#         rmin=np.min(r[ilist==i])
#         rmax=np.max(r[ilist==i])
#         rmed_mr=np.append(rmed_mr,np.median(r[ilist==i]))
#         umed_mr=np.append(umed_mr,np.median(gu[ikeep][ilist==i]))
#         rhomed_mr=np.append(rhomed_mr,np.median(grho[ikeep][ilist==i]))
#         rhoshell_mr=np.append(rhoshell_mr,3*np.sum(gmass[ikeep][ilist==i])/4./np.pi/(rmax**3-rmin**3))

        
        
