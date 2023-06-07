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
    '''
    Read simulation snapshots - can be in HDF5 or GADGET binary format.
    
        Parameters:
            snapfilename - path of the snapshot
            snapfileformat - format of the snapshot; can be 1 or 2 (binary)
            or 3 (HDF5), or 'SNAP1/SNAP2' (binary) of 'HDF5'
            convention - 'SWIFT', 'GADGET4', otherwise assumes 'GADGET2/3'
        
        If the file is in HDF5 formay, omit the .hdf5 suffix.
    '''
    def __init__(self,snapfilename,snapfileformat,**kwargs):
        self.snapfilename=snapfilename
        self.snapfileformat=snapfileformat
        if self.snapfileformat=='HDF5' or self.snapfileformat=='3':
            self.convention=kwargs.get('convention')
        
    def test(self):
        filename=self.snapfilename+'.hdf5'
        if os.path.exists(filename)==False:
            filename=self.snapfilename+'.0.hdf5'
        print('Reading data from %s'%filename)
        
        with h5py.File(filename,'r') as f:
            self.NumPart_Total=f['Header'].attrs['NumPart_Total'][()]
            # part=ParticleProperties(1,self.NumPart_Total[1])

        
    def ReadSnapshot(self):
        '''
        Reads in data from a single or multiple snapshots.
        '''

        self.NumPartType=6
        
        if self.snapfileformat=='HDF5' or self.snapfileformat=='3':
            filename=self.snapfilename+'.hdf5'
            if os.path.exists(filename)==False:
              filename=self.snapfilename+'.0.hdf5'
            print('Reading data from %s'%filename)
            
            with h5py.File(filename,'r') as f:
                self.ScaleFactor=f['Header'].attrs['Time']
                NumFiles=f['Header'].attrs['NumFilesPerSnapshot']
                self.NumPart_Total=f['Header'].attrs['NumPart_Total'][()]
                self.NumPartType=np.shape(self.NumPart_Total)[0]
                if self.convention=='SWIFT':
                    self.BoxSize=f['Header'].attrs['BoxSize'][()][0]
                    self.OmegaDM=f['Cosmology'].attrs['Omega_cdm']
                    self.OmegaBar=f['Cosmology'].attrs['Omega_b']
                    self.Omega0=self.OmegaDM+self.OmegaBar
                    self.OmegaLambda=f['Cosmology'].attrs['Omega_lambda']
                    self.HubbleParam=f['Cosmology'].attrs['h']
                elif self.convention=='GADGET4':
                    self.BoxSize=f['Header'].attrs['BoxSize']
                    self.Omega0=f['Parameters'].attrs['Omega0']
                    self.OmegaLambda=f['Parameters'].attrs['OmegaLambda']
                    self.HubbleParam=f['Parameters'].attrs['HubbleParam']
                else:
                    self.BoxSize=f['Header'].attrs['BoxSize']
                    self.Omega0=f['Header'].attrs['Omega0']
                    self.OmegaLambda=f['Header'].attrs['OmegaLambda']
                    self.HubbleParam=f['Header'].attrs['HubbleParam']
                if NumFiles>1:
                    print('Data is split across %d files'%NumFiles)
                
            NumPart=np.sum(self.NumPart_Total)
            print('Number of particles: %010d'%NumPart)
            print('Number of particle types: %d'%self.NumPartType)
            self.pos=np.ndarray(shape=(NumPart,3))
            self.vel=np.ndarray(shape=(NumPart,3))
            self.pids=np.ndarray(shape=(NumPart),dtype=np.uint64)
            if self.NumPart_Total[0]>0:
                self.u=np.ndarray(shape=(NumPart_Total[0]))
                self.rho=np.ndarray(shape=(NumPart_Total[0]))
                
            istart=np.zeros(self.NumPartType,dtype=np.uint64)
            offset=0
            
            for i in range(1,self.NumPartType):
                offset+=self.NumPart_Total[i-1]
                if self.NumPart_Total[i]>0:
                    istart[i]=offset
            ifinish=np.copy(istart)
            
            if NumFiles>1:
                for i in range(NumFiles):
                    filename=self.snapfilename+'.%d.hdf5'%i
                    print('Reading in file %s...'%filename)
                    with h5py.File(filename,'r') as f:
                        NumPart_ThisFile=f['Header'].attrs['NumPart_ThisFile'][()]

                        for itype in range(self.NumPartType):
                            if NumPart_ThisFile[itype]>0:
                                ifinish[itype]=istart[itype]+NumPart_ThisFile[itype]
                                self.pos[istart[itype]:ifinish[itype]]=f['PartType%d/Coordinates'%itype][()]
                                self.vel[istart[itype]:ifinish[itype]]=f['PartType%d/Velocities'%itype][()]
                                self.pids[istart[itype]:ifinish[itype]]=f['PartType%d/ParticleIDs'%itype][()]
                                if self.NumPart_Total[0]>0:
                                    self.u[istart[itype]:ifinish[itype]]=f['PartType%d/InternalEnergy'%itype][()]
                                    self.rho[istart[itype]:ifinish[itype]]=f['PartType%d/Density'%itype][()]
                                istart[itype]=ifinish[itype]
            else:
                with h5py.File(filename,'r') as f:
                    for itype in range(self.NumPartType):
                        if self.NumPart_Total[itype]>0:
                            ifinish[itype]=istart[itype]+self.NumPart_Total[itype]
                            self.pos[istart[itype]:ifinish[itype]]=f['PartType%d/Coordinates'%itype][()]
                            self.vel[istart[itype]:ifinish[itype]]=f['PartType%d/Velocities'%itype][()]
                            self.pids[istart[itype]:ifinish[itype]]=f['PartType%d/ParticleIDs'%itype][()]
                            if self.NumPart_Total[0]>0:
                                self.u[istart[itype]:ifinish[itype]]=f['PartType%d/InternalEnergy'%itype][()]
                                self.rho[istart[itype]:ifinish[itype]]=f['PartType%d/Density'%itype][()]
                            istart[itype]=ifinish[itype]
        else:
            fileroot=self.snapfilename
            
            filename=fileroot
            
            if os.path.exists(filename)==False:
              filename=fileroot+'.0'
            print('Reading data from %s'%filename)

            with open(filename,'rb') as f:
                offset=0

                if self.snapfileformat=='SNAP2':
                    offset+=16

                offset+=4
                f.seek(offset,os.SEEK_SET)
                self.NumPart_ThisFile=np.fromfile(f,dtype=np.int32,count=6)

                offset+=24
                f.seek(offset,os.SEEK_SET)
                self.MassTable=np.fromfile(f,dtype=np.float64,count=6)

                offset+=48
                f.seek(offset,os.SEEK_SET)
                self.ScaleFactor=np.fromfile(f,dtype=np.float64,count=1)[0]

                offset+=8
                f.seek(offset,os.SEEK_SET)
                self.Redshift=np.fromfile(f,dtype=np.float64,count=1)[0]

                offset+=16
                f.seek(offset,os.SEEK_SET)
                self.NumPart_Total=np.fromfile(f,dtype=np.int32,count=6)
    
                offset+=28
                f.seek(offset,os.SEEK_SET)
                NumFiles=np.fromfile(f,dtype=np.int32,count=1)[0]

                offset+=4
                f.seek(offset,os.SEEK_SET)
                self.BoxSize=np.fromfile(f,dtype=np.float64,count=1)[0]

                offset+=8
                f.seek(offset,os.SEEK_SET)
                self.Omega0=np.fromfile(f,dtype=np.float64,count=1)[0]

                offset+=8
                f.seek(offset,os.SEEK_SET)
                self.OmegaLambda=np.fromfile(f,dtype=np.float64,count=1)[0]

                offset+=8
                f.seek(offset,os.SEEK_SET)
                self.HubbleParam=np.fromfile(f,dtype=np.float64,count=1)[0]

                if NumFiles>1:
                    print('Data is split across %d files'%NumFiles)
            f.close()

            NumPart=np.sum(self.NumPart_Total)
            print('Number of particles: %010d'%NumPart)
            self.pos=np.ndarray(shape=(NumPart,3))
            self.vel=np.ndarray(shape=(NumPart,3))
            self.pids=np.ndarray(shape=(NumPart),dtype=np.uint64)
            if self.NumPart_Total[0]>0:
                self.u=np.ndarray(shape=(self.NumPart_Total[0]))
                self.rho=np.ndarray(shape=(self.NumPart_Total[0]))
                
            istart=np.zeros(self.NumPartType,dtype=np.uint64)
            offset=0

            for i in range(1,self.NumPartType):
                offset+=self.NumPart_Total[i-1]
                if self.NumPart_Total[i]>0:
                    istart[i]=offset
            ifinish=np.copy(istart)

            for ifile in range(NumFiles):
                if ifile>0:
                    filename=fileroot+'.%d'%ifile
                with open(filename,'rb') as f:
                    offset=0
                    if self.snapfileformat=='SNAP2':
                        offset+=16
                    offset+=4
                    f.seek(offset,os.SEEK_SET)
                    NumPartInThisFile=np.fromfile(f,dtype=np.int32,count=6)

                    NumPartInFile=np.sum(NumPartInThisFile)

                    print('Reading %010d particles from %s'%(NumPartInFile,filename))

                    offset=264   # includes 2 x 4 byte buffers

                    if self.snapfileformat=='SNAP2':
                        offset+=16
        
                    offset+=4    # 1st 4 byte buffer
                    if self.snapfileformat=='SNAP2':
                        offset+=16

                    f.seek(offset,os.SEEK_SET)
                    
                    pos_block=np.fromfile(f,dtype=np.float32,count=3*NumPartInFile)

                    # Increment beyond the POS block
                    offset+=3*NumPartInFile*4
                    offset+=4   # 2nd 4 byte buffer

                    # Open the VEL block
                    offset+=4   # 1st 4 byte buffer
                    if self.snapfileformat=='SNAP2':
                        offset+=16

                    f.seek(offset,os.SEEK_SET)
                    
                    vel_block=np.fromfile(f,dtype=np.float32,count=3*NumPartInFile)

                    # Increment beyond the VEL block
                    offset+=3*NumPartInFile*4
                    offset+=4   # 2nd 4 byte buffer

                    # Open the IDS block
                    offset+=4   # 1st 4 byte buffer
                    if self.snapfileformat=='SNAP2':
                        offset+=16

                    f.seek(offset,os.SEEK_SET)

                    pids_block=np.fromfile(f,dtype=np.uint64,count=NumPartInFile)
    
#                    print(pids_block)

#                    # Copy these blocks into pos,vel,pids arrays
#                    SingleOffsetStart=np.zeros(self.NumPartType,dtype=np.uint64)
#                    TripleOffsetStart=np.zeros(self.NumPartType,dtype=np.uint64)
#                    TypeOffset=0
#                    for i in range(1,self.NumPartType):
#                        TypeOffset+=self.NumPart_ThisFile[i-1]
#                        if self.NumPart_ThisFile[i]>0:
#                            SingleOffsetStart[i]=TypeOffset
#                            TripleOffsetStart[i]=TypeOffset
#                    SingleOffsetFinish=np.copy(SingleOffsetStart)
#                    TripleOffsetFinish=np.copy(TripleOffsetStart)
                    ifinish=istart+NumPartInThisFile.astype(np.uint64)
                    
                    astart=0
                    bstart=0
                    
                    for itype in range(self.NumPartType):
                        if NumPartInThisFile[itype]>0:
                            afinish=astart+3*NumPartInThisFile[itype]
                            bfinish=bstart+NumPartInThisFile[itype]
                            self.pos[istart[itype]:ifinish[itype]]=pos_block[astart:afinish].reshape(NumPartInThisFile[itype],3)
                            self.vel[istart[itype]:ifinish[itype]]=vel_block[astart:afinish].reshape(NumPartInThisFile[itype],3)
                            print(pids_block.size)
                            self.pids[istart[itype]:ifinish[itype]]=pids_block[bstart:bfinish]
                            astart=afinish
                            bstart=bfinish
#                            SingleOffsetFinish[itype]=SingleOffsetStart[itype]+self.NumPart_ThisFile[itype]
#                            TripleOffsetFinish[itype]=TripleOffsetStart[itype]+3*self.NumPart_ThisFile[itype]
#                            self.pos[istart[itype]:ifinish[itype]]=pos_block[TripleOffsetStart[itype]:TripleOffsetFinish[itype]].reshape(self.NumPart_ThisFile[itype],3)
#                            self.vel[istart[itype]:ifinish[itype]]=vel_block[TripleOffsetStart[itype]:TripleOffsetFinish[itype]].reshape(self.NumPart_ThisFile[itype],3)
#                            self.pids[istart[itype]:ifinish[itype]]=pids_block[SingleOffsetStart[itype]:SingleOffsetFinish[itype]]
                    istart=ifinish
#                            SingleOffsetStart[itype]=SingleOffsetFinish[itype]
#                            TripleOffsetStart[itype]=TripleOffsetFinish[itype]

                    print('ifinish',ifinish)
#                    print('NumPartInThisFile',NumPartInThisFile)
#                    # Increment beyond the IDS block
#                    offset+=NumPart*8
#                    offset+=4   # 2nd 4 byte buffer
#
#                    offset+=4   # 1st 4 byte buffer
#                    if self.snapfileformat=='SNAP2':
#                    offset+=16
#
#                    f.seek(offset,os.SEEK_SET)
#
#                    TypeOffsetStart=np.zeros(self.NumPartType,dtype=np.uint64)
#                    TypeOffset=0
#                    for i in range(1,self.NumPartType):
#                        TypeOffset+=self.NumPart_ThisFile[i-1]
#                        if self.NumPart_ThisFile[i]>0:
#                            TypeOffsetStart[i]=TypeOffset
#                    print(TypeOffsetStart)
#                    print(self.NumPart_ThisFile)
#                    TypeOffsetFinish=np.copy(TypeOffsetStart)
#
#                    pids_block=np.fromfile(f,dtype=np.uint64,count=NumPart)
#
#                    for itype in range(self.NumPartType):
#                        if self.NumPart_ThisFile[itype]>0:
#                            ifinish[itype]=istart[itype]+self.NumPart_ThisFile[itype]
#                            TypeOffsetFinish[itype]=TypeOffsetStart[itype]+self.NumPart_ThisFile[itype]
#                            self.pids[istart[itype]:ifinish[itype]]=pids_block[TypeOffsetStart[itype]:TypeOffsetFinish[itype]].reshape(self.NumPart_ThisFile[itype],3)
#                            istart[itype]=ifinish[itype]
#                            TypeOffsetStart[itype]=TypeOffsetFinish[itype]
#


#    def GetParticlesByType
#    def GetParticleSubset(self,subset_pids):
#        elements_in_particle_ids=np.where(np.in1d(self.pids,subset_pids))[0]
#        return elements_in_particle_ids
