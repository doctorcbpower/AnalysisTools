#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 13:55:13 2023

@author: cpower

Script to read/write various snapshot data

"""
import h5py
import numpy as np
import os

class SnapshotTools:
    def __init__(self,snapfilename,snapfileformat):
        self.snapfilename=snapfilename
        self.snapfileformat=snapfileformat
            
    def ReadSnapshot(self):        
        if self.snapfileformat=='HDF5':
            filename=self.snapfilename+'.hdf5'
            if os.path.exists(filename)==False:
              filename=self.snapfilename+'.0.hdf5'            
            print('Reading data from %s'%filename)
            with h5py.File(filename,'r') as f:
                self.ScaleFactor=f['Header'].attrs['Time']
                self.BoxSize=f['Header'].attrs['BoxSize']
                NumFiles=f['Header'].attrs['NumFilesPerSnapshot']
                self.NumPart_Total=f['Header'].attrs['NumPart_Total'][()]
                if NumFiles>1:
                    print('Data is split across %d files'%NumFiles)
                
            self.pos=np.ndarray(shape=(self.NumPart_Total[1],3))
            self.vel=np.ndarray(shape=(self.NumPart_Total[1],3))
            self.pids=np.ndarray(shape=(self.NumPart_Total[1]))

            if NumFiles>1:
                istart=0
                for i in range(NumFiles):                
                    filename=self.snapfilename+'.%d.hdf5'%i
                    print('Reading in file %s...'%filename)
                    with h5py.File(filename,'r') as f:
                        NumPart_ThisFile=f['Header'].attrs['NumPart_ThisFile'][()]                        
                        ifinish=istart+NumPart_ThisFile[1]                        
                        self.pos[istart:ifinish]=f['PartType1/Coordinates'][()]
                        self.vel[istart:ifinish]=f['PartType1/Velocities'][()]
                        self.pids[istart:ifinish]=f['PartType1/ParticleIDs'][()]                        
                        istart=ifinish
            else:                
                with h5py.File(filename,'r') as f:
                    self.pos=f['PartType1/Coordinates'][()]
                    self.vel=f['PartType1/Coordinates'][()]
                    self.pids=f['PartType1/ParticleIDs'][()]
                
    def GetParticleSubset(self,subset_pids):
        elements_in_particle_ids=np.where(np.in1d(self.pids,subset_pids))[0]
        return elements_in_particle_ids
    
