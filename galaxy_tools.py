#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 13:55:13 2023

@author: cpower

Script to read and analyse halo catalogues generated using SubFind and VELOCIraptor.
"""
import h5py
import numpy as np
import os

class GalaxyTools:
    '''
    Read SHARK catalogues
    
        Parameters:
            galfilename - path of the galaxy catalogue
    '''
    def __init__(self,galfilename,galfileformat,**kwargs):
        self.galfilename=galfilename
        self.galfileformat=galfileformat

    def ReadGalaxyCatalogue(self):
        '''
        Reads in data from a single or multiple snapshots.
        '''
        
        print('Reading data from %s'%self.galfilename)
            
        with h5py.File(self.galfilename,'r') as f:
            print(f['galaxies'].keys())
            self.pos=np.array([f['galaxies']['position_x'][()],\
                               f['galaxies']['position_y'][()],\
                               f['galaxies']['position_z'][()]]).T
            self.vel=np.array([f['galaxies']['velocity_x'][()],\
                               f['galaxies']['velocity_y'][()],\
                               f['galaxies']['velocity_z'][()]]).T
            self.vel=np.array([f['galaxies']['velocity_x'][()],\
                               f['galaxies']['velocity_y'][()],\
                               f['galaxies']['velocity_z'][()]]).T
            self.mstar_disk=f['galaxies']['mstars_disk'][()]
            self.mstar_bulge=f['galaxies']['mstars_bulge'][()]
            self.mgas_disk=f['galaxies']['mgas_disk'][()]
            self.mgas_bulge=f['galaxies']['mgas_bulge'][()]
            self.matom_disk=f['galaxies']['matom_disk'][()]
            self.matom_bulge=f['galaxies']['matom_bulge'][()]
            self.mstar_tot=self.mstar_disk+self.mstar_bulge
            self.mgas_tot=self.mgas_disk+self.mgas_bulge
            self.matom_tot=self.matom_disk+self.matom_bulge

