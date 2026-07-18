#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:45:15 2023

@author: 00075868
"""

import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
from . import snapshot_tools as st

class ProfileTools:
    """
    General purpose particle profile analysis.

    Usage: 
    from analysis_tools import ProfileTools

    pt = ProfileTools(numbins=40)

    rho = pt.volume_density(
        pos,
        mass,
        centre,
        rmin=0.1,
        rmax=300
    )
    
    sigma = pt.surface_density(
            pos,
            mass,
            centre,
            rmin=0.1,
            rmax=50
    )

    hz = pt.scale_height(
            pos,
            mass,
            centre,
            rmin=0.1,
            rmax=50
    )

    pt.plot(
        rho,
        "density",
        ylabel=r"$\rho(r)$"
    )
    """

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
            
    ## Compute positions
    @staticmethod
    def relative_position(pos, centre=None):

        if centre is None:
            centre = np.zeros(3)

        return pos - centre


    def spherical_radius(self, pos, centre=None):

        dpos = self.relative_position(pos, centre)

        return np.sqrt(np.sum(dpos**2, axis=1))


    def cylindrical_radius(self, pos, centre=None):

        dpos = self.relative_position(pos, centre)

        return np.sqrt(dpos[:,0]**2 + dpos[:,1]**2)


    ## Impose binning
    def radial_bins(self, rmin, rmax, nbins=None,
                    logarithmic=True):

        if nbins is None:
            nbins = self.numbins

        if logarithmic:
            return np.logspace(
                np.log10(rmin),
                np.log10(rmax),
                nbins+1
            )

        return np.linspace(rmin, rmax, nbins+1)


    def bin_indices(self, r, bins):
        return np.digitize(r, bins)-1


    ## Calculate radial profiles
    def profile(self, r, quantity, bins,
                weights=None):

        index = self.bin_indices(r, bins)

        nbins = len(bins)-1

        result = {
            "r": np.zeros(nbins),
            "mean": np.full(nbins,np.nan),
            "median": np.full(nbins,np.nan),
            "count": np.zeros(nbins)
        }


        for i in range(nbins):

            mask = index == i

            result["r"][i] = 0.5*(bins[i]+bins[i+1])

            if np.any(mask):

                q = quantity[mask]

                if weights is None:
                    result["mean"][i] = np.mean(q)
                else:
                    result["mean"][i] = np.average(
                        q,
                        weights=weights[mask]
                    )

                result["median"][i] = np.median(q)
                result["count"][i] = np.sum(mask)

        return result

    ## Calculate density profiles
    def volume_density(self, pos, mass, centre, rmin, rmax, nbins=None):
        r = self.spherical_radius(pos, centre)

        bins = self.radial_bins(rmin,rmax,nbins)

        index = self.bin_indices(r,bins)

        rho=[]
        radius=[]

        for i in range(len(bins)-1):

            mask=index==i

            volume = (4*np.pi/3.) * ( bins[i+1]**3-bins[i]**3)

            rho.append(np.sum(mass[mask])/volume)

            radius.append(0.5*(bins[i]+bins[i+1]))

        return {
            "r":np.array(radius),
            "density":np.array(rho)
        }

    def surface_density(self,pos,mass, centre, rmin,rmax, nbins=None):
        R=self.cylindrical_radius(pos,centre)

        bins=self.radial_bins(rmin,rmax,nbins)

        index=self.bin_indices(R,bins)

        sigma=[]
        radius=[]

        for i in range(len(bins)-1):

            mask=index==i

            area=np.pi*(bins[i+1]**2-bins[i]**2)

            sigma.append(np.sum(mass[mask])/area)

            radius.append(0.5*(bins[i]+bins[i+1]))

        return {
            "r":np.array(radius),
            "density":np.array(sigma)
        }


    ## Calculate kinematic profiles
    def velocity_dispersion(self,pos,vel, centre, rmin,rmax, nbins=None):
        r=self.spherical_radius(pos,centre)

        bins=self.radial_bins(rmin,rmax,nbins)

        index=self.bin_indices(r,bins)

        sigma=[]

        radius=[]

        for i in range(len(bins)-1):

            mask=index==i

            radius.append(0.5*(bins[i]+bins[i+1]))

            if np.sum(mask)>5:
                v=vel[mask]
                mean=np.mean(v,axis=0)
                sigma.append(np.sqrt(np.mean(np.sum(v*v,axis=1))-np.sum(mean*mean)))
            else:
                sigma.append(np.nan)

        return {
            "r":np.array(radius),
            "sigma":np.array(sigma)
        }


    ## Calculate disc structure
    def scale_height(self,pos,mass,
                     centre,
                     rmin,rmax,
                     nbins=None):

        R=self.cylindrical_radius(pos,centre)

        bins=self.radial_bins(rmin,rmax,nbins)

        index=self.bin_indices(R,bins)

        hz=[]
        radius=[]

        z=pos[:,2]-centre[2]

        for i in range(len(bins)-1):

            mask=index==i

            radius.append(0.5*(bins[i]+bins[i+1]))

            if np.any(mask):
                hz.append(np.sqrt(np.average(z[mask]**2,weights=mass[mask])))
            else:
                hz.append(np.nan)

        return {
            "r":np.array(radius),
            "scale_height":np.array(hz)
        }

    def vertical_velocity_dispersion(self,pos,vel, centre, rmin,rmax, nbins=None):
        R=self.cylindrical_radius(pos,centre)

        bins=self.radial_bins(rmin,rmax,nbins)

        index=self.bin_indices(R,bins)

        sigma_z=[]

        radius=[]
        
        vz = vel[:,2]

        for i in range(len(bins)-1):

            mask=index==i

            radius.append(0.5*(bins[i]+bins[i+1]))

            if np.sum(mask)>5:
                v=vz[mask]
                mean=np.mean(v,axis=0)
                sigma_z.append(np.sqrt(np.mean(v*v)-np.sum(mean*mean)))
            else:
                sigma_z.append(np.nan)

        return {
            "r":np.array(radius),
            "sigma_z":np.array(sigma_z)
        }


    # ------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------

    def plot(self, profile,ykey,xlabel="Radius",ylabel=None,xlog=True,ylog=True,label=None):
        plt.figure()

        plt.plot(
            profile["r"],
            profile[ykey],
            label=label
        )

        if xlog:
            plt.xscale("log")
        if ylog:
            plt.yscale("log")

        plt.xlabel(xlabel)

        if ylabel:
            plt.ylabel(ylabel)

        if label:
            plt.legend()

        plt.show()


class MassFunctionTools:
    def __init__(self,**kwargs):
        self.comoving_units=True
        # Use comoving or physical units when plotting?
        if 'comoving_units' in kwargs:
            self.comoving_units=kwargs.get('comoving_units')
        self.numbins=25
        # Number of bins to plot
        if 'numbins' in kwargs:
            self.numbins=kwargs.get('numbins')
        self.halo_id=0

    def BinByHaloMass(self,
                      halo_mass,
                      lbox,
                      delta_logmass=0.,
                      numbins=10,
                      centre_stat="median"):     # "median" or "mean"):

        logm = np.log10(halo_mass)

        # Case 1: user specifies delta_logmass
        if delta_logmass > 0:
            dm = logm.max() - logm.min()
            numbins = max(1, int(np.ceil(dm / delta_logmass)))
            dlm = delta_logmass
        else:
            # Case 2: user specifies numbins
            dm = logm.max() - logm.min()
            dlm = dm / numbins

        print(f"Binning data with {numbins} bins, width dlm = {dlm:.4f} dex")

        # -----------------------------------------
        # Define bin edges from dlm
        # -----------------------------------------
        bin_edges = np.arange(logm.min(), logm.max() + dlm, dlm)
        numbins = len(bin_edges) - 1  # update

        # Digitize each mass
        bin_index = np.digitize(logm, bin_edges) - 1

        # Remove values outside bin range
        valid = (bin_index >= 0) & (bin_index < numbins)
        bin_index = bin_index[valid]
        logm_valid = logm[valid]

        # Count halos in each bin
        num = np.bincount(bin_index, minlength=numbins)

        # Avoid zero counts
        num_safe = np.where(num > 0, num, np.nan)

        # -----------------------------------------
        # Compute lm as mean or median mass in each bin
        # -----------------------------------------
        lm = np.full(numbins, np.nan)
        for i in range(numbins):
            mask = bin_index == i
            if np.any(mask):
                if centre_stat.lower() == "mean":
                    lm[i] = np.mean(logm_valid[mask])
                else:
                    lm[i] = np.median(logm_valid[mask])  # default

        # dn/dlogM (log10)
        ldndlm = np.log10(num_safe / lbox**3 / dlm)

        return lm, ldndlm
