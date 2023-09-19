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
    def __init__(self,snapfilename,snapfileformat,gas_type=0,dm_type=1,star_type=4,bh_type=5,**kwargs):
        self.snapfilename=snapfilename
        self.snapfileformat=snapfileformat
        self.gas_type=gas_type
        self.dm_type=dm_type
        self.star_type=star_type
        self.bh_type=bh_type
        if self.snapfileformat=='HDF5' or self.snapfileformat=='3':
            self.convention=kwargs.get('convention')
        
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
                NumFiles=f['Header'].attrs['NumFilesPerSnapshot']
                NameOfMassBlock='Mass'
                if 'Masses' in list(f['PartType1'].keys()):
                    NameOfMassBlock='Masses'
                self.NumPart_Total=f['Header'].attrs['NumPart_Total'][()]
                self.NumPartType=np.shape(self.NumPart_Total)[0]
                if self.convention=='SWIFT':
                    self.ScaleFactor=f['Header'].attrs['Scale-factor']
                    self.BoxSize=f['Header'].attrs['BoxSize'][()][0]
                    self.OmegaDM=f['Cosmology'].attrs['Omega_cdm']
                    self.OmegaBar=f['Cosmology'].attrs['Omega_b']
                    self.Omega0=self.OmegaDM+self.OmegaBar
                    self.OmegaLambda=f['Cosmology'].attrs['Omega_lambda']
                    self.HubbleParam=f['Cosmology'].attrs['h']
                    self.MassTable=f['Header'].attrs['MassTable'][()]
                elif self.convention=='GADGET4':
                    self.ScaleFactor=f['Header'].attrs['Time']
                    self.BoxSize=f['Header'].attrs['BoxSize']
                    self.Omega0=f['Parameters'].attrs['Omega0']
                    self.OmegaLambda=f['Parameters'].attrs['OmegaLambda']
                    self.HubbleParam=f['Parameters'].attrs['HubbleParam']
                    self.MassTable=f['Header'].attrs['MassTable'][()]
                else:
                    self.ScaleFactor=f['Header'].attrs['Time']
                    self.BoxSize=f['Header'].attrs['BoxSize']
                    self.Omega0=f['Header'].attrs['Omega0']
                    self.OmegaLambda=f['Header'].attrs['OmegaLambda']
                    self.HubbleParam=f['Header'].attrs['HubbleParam']
                    self.MassTable=f['Header'].attrs['MassTable'][()]
                print('Simulation scale factor: %lf'%self.ScaleFactor)
                if NumFiles>1:
                    print('Data is split across %d files'%NumFiles)
                
            NumPart=np.sum(self.NumPart_Total)
            print('Number of particles: %010d'%NumPart)
            print('Number of particle types: %d'%self.NumPartType)
            idx_with_mass=np.where(self.MassTable>0)[0]
            NumPart_InMassBlock=np.sum(self.NumPart_Total[idx_with_mass])
            if NumPart_InMassBlock>0:
                print('Number of particles in mass block: %010d'%NumPart_InMassBlock)

            self.pos=np.ndarray(shape=(NumPart,3))
            self.vel=np.ndarray(shape=(NumPart,3))
            self.pids=np.ndarray(shape=(NumPart),dtype=np.uint64)
            self.mass=np.ndarray(shape=(NumPart))
            
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
                                if self.MassTable[itype]==0:
                                    if NameOfMassBlock=='Mass':
                                        self.mass[istart[itype]:ifinish[itype]]=f['PartType%d/Mass'%itype][()]
                                    else:
                                        if itype==5 and self.convention=='SWIFT':
                                            self.mass[istart[itype]:ifinish[itype]]=f['PartType%d/DynamicalMasses'%itype][()]
                                        else:
                                           self.mass[istart[itype]:ifinish[itype]]=f['PartType%d/Masses'%itype][()]
                                if itype==0:
                                    if self.convention=='SWIFT':
                                        self.u[istart[itype]:ifinish[itype]]=f['PartType%d/InternalEnergies'%itype][()]
                                        self.rho[istart[itype]:ifinish[itype]]=f['PartType%d/Densities'%itype][()]
                                    else:
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
                            if self.MassTable[itype]==0:
                                if NameOfMassBlock=='Mass':
                                    self.mass[istart[itype]:ifinish[itype]]=f['PartType%d/Mass'%itype][()]
                                else:
                                    if itype==5 and self.convention=='SWIFT':
                                        self.mass[istart[itype]:ifinish[itype]]=f['PartType%d/DynamicalMasses'%itype][()]
                                    else:
                                       self.mass[istart[itype]:ifinish[itype]]=f['PartType%d/Masses'%itype][()]
                            if itype==0:
                                if self.convention=='SWIFT':
                                    self.u[istart[itype]:ifinish[itype]]=f['PartType%d/InternalEnergies'%itype][()]
                                    self.rho[istart[itype]:ifinish[itype]]=f['PartType%d/Densities'%itype][()]
                                else:
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

            idx_with_mass=np.where(self.MassTable==0)[0]    # Want to know which species are in the mass block, so
                                                            # their MassTable entries will be zero
            NumPart=np.sum(self.NumPart_Total)
            print('Number of particles: %010d'%NumPart)
            NumPart_InMassBlock=np.sum(self.NumPart_Total[idx_with_mass])
            print('Number of particles in mass block: %010d'%NumPart_InMassBlock)
            
            self.pos=np.ndarray(shape=(NumPart,3))
            self.vel=np.ndarray(shape=(NumPart,3))
            self.pids=np.ndarray(shape=(NumPart),dtype=np.uint64)
            self.mass=np.ndarray(shape=(NumPart))
            
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

                    NumPart_InMassBlock_InFile=np.sum(NumPartInThisFile[idx_with_mass])
                    
                    if NumFiles>1:
                        print('Reading %010d particles from %s'%(NumPartInFile,filename))
                        print('Number of particles in mass block: %010d'%NumPart_InMassBlock_InFile)

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
                    
                   # Increment beyond the IDs block
                    offset+=NumPartInFile*4
                    offset+=4   # 2nd 4 byte buffer

                    # Open the mass block
                    offset+=4   # 1st 4 byte buffer
                    if self.snapfileformat=='SNAP2':
                        offset+=16

                    f.seek(offset,os.SEEK_SET)

                    mass_block=np.fromfile(f,dtype=np.float32,count=NumPart_InMassBlock_InFile)

                    if NumPartInThisFile[0]>0:
                        # Increment beyond the mass block
                        offset+=NumPart_InMassBlock_InFile*4
                        offset+=4   # 2nd 4 byte buffer

                        # Open the mass block
                        offset+=4   # 1st 4 byte buffer
                        if self.snapfileformat=='SNAP2':
                            offset+=16

                        f.seek(offset,os.SEEK_SET)

                        u_block=np.fromfile(f,dtype=np.float32,count=NumPartInThisFile[0])

                       # Increment beyond the mass block
                        offset+=NumPartInThisFile[0]*4
                        offset+=4   # 2nd 4 byte buffer

                        # Open the mass block
                        offset+=4   # 1st 4 byte buffer
                        if self.snapfileformat=='SNAP2':
                            offset+=16

                        f.seek(offset,os.SEEK_SET)

                        rho_block=np.fromfile(f,dtype=np.float32,count=NumPartInThisFile[0])


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
                    cstart=0
                    
                    for itype in range(self.NumPartType):
                        if NumPartInThisFile[itype]>0:
                            afinish=astart+3*NumPartInThisFile[itype]
                            bfinish=bstart+NumPartInThisFile[itype]
                            self.pos[istart[itype]:ifinish[itype]]=pos_block[astart:afinish].reshape(NumPartInThisFile[itype],3)
                            self.vel[istart[itype]:ifinish[itype]]=vel_block[astart:afinish].reshape(NumPartInThisFile[itype],3)
                            self.pids[istart[itype]:ifinish[itype]]=pids_block[bstart:bfinish]
                            if np.isin(itype,idx_with_mass)==True:
                                cfinish=cstart+NumPartInThisFile[itype]
                                self.mass[istart[itype]:ifinish[itype]]=mass_block[cstart:cfinish]
                            else:
                                self.mass[istart[itype]:ifinish[itype]]=self.MassTable[itype]*np.ones(NumPartInThisFile[itype])

                            self.u[istart[itype]:ifinish[itype]]=u_block[bstart:bfinish]
                            self.rho[istart[itype]:ifinish[itype]]=rho_block[bstart:bfinish]
                            astart=afinish
                            bstart=bfinish
                            cstart=cfinish
#                            SingleOffsetFinish[itype]=SingleOffsetStart[itype]+self.NumPart_ThisFile[itype]
#                            TripleOffsetFinish[itype]=TripleOffsetStart[itype]+3*self.NumPart_ThisFile[itype]
#                            self.pos[istart[itype]:ifinish[itype]]=pos_block[TripleOffsetStart[itype]:TripleOffsetFinish[itype]].reshape(self.NumPart_ThisFile[itype],3)
#                            self.vel[istart[itype]:ifinish[itype]]=vel_block[TripleOffsetStart[itype]:TripleOffsetFinish[itype]].reshape(self.NumPart_ThisFile[itype],3)
#                            self.pids[istart[itype]:ifinish[itype]]=pids_block[SingleOffsetStart[itype]:SingleOffsetFinish[itype]]
                    istart=ifinish
#                            SingleOffsetStart[itype]=SingleOffsetFinish[itype]
#                            TripleOffsetStart[itype]=TripleOffsetFinish[itype]

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
    def LoadParticlesByType(self,part_type='all'):
        isgas=np.sum(self.NumPart_Total[:self.gas_type])
        ifgas=np.sum(self.NumPart_Total[:self.gas_type+1])
        isdm=np.sum(self.NumPart_Total[:self.dm_type])
        ifdm=np.sum(self.NumPart_Total[:self.dm_type+1])
        isstar=np.sum(self.NumPart_Total[:self.star_type])
        ifstar=np.sum(self.NumPart_Total[:self.star_type+1])
        isbh=np.sum(self.NumPart_Total[:self.bh_type])
        ifbh=np.sum(self.NumPart_Total[:self.bh_type+1])

        loadgas=False
        loaddm=True
        loadstar=False
        loadbh=False

        if part_type=='all':
            loadgas=True
            loadstar=True
            loadbh=True
        if part_type=='gas':
            loadgas=True
            loaddm=False
        if part_type=='star':
            loadstar=True
            loaddm=False
        if part_type=='bh':
            loadbh=True
            loaddm=False
        
        if loadgas==True:
            self.gas=self.ParticleProperties(self.pos[isgas:ifgas],
                                        self.vel[isgas:ifgas],
                                        self.pids[isgas:ifgas],
                                        self.mass[isgas:ifgas],
                                        internal_energy=self.u[isgas:ifgas],
                                        density=self.rho[isgas:ifgas])

        if loaddm==True:
            self.dm=self.ParticleProperties(self.pos[isdm:ifdm],
                                        self.vel[isdm:ifdm],
                                        self.pids[isdm:ifdm],
                                        self.mass[isdm:ifdm])

        if loadstar==True:
            self.star=self.ParticleProperties(self.pos[isstar:ifstar],
                                        self.vel[isstar:ifstar],
                                        self.pids[isstar:ifstar],
                                        self.mass[isstar:ifstar])
        
        if loadbh==True:
            self.bh=self.ParticleProperties(self.pos[isbh:ifbh],
                                        self.vel[isbh:ifbh],
                                        self.pids[isbh:ifbh],
                                        self.mass[isbh:ifbh])
        
    class ParticleProperties:
        def __init__(self,pos,vel,pids,mass,**kwargs):
            self.pos=pos
            self.vel=vel
            self.pids=pids
            self.mass=mass
            if 'internal_energy' in kwargs:
                self.internal_energy=kwargs.get('internal_energy')
            if 'density' in kwargs:
                self.density=kwargs.get('density')
                

    def UnitConversion(self,**kwargs):
        if kwargs.get('convert_to_physical')!=None:
            self.convert_to_physical=kwargs.get('convert_to_physical')
            if self.convert_to_physical==True:
                self.pos*=self.ScaleFactor
                self.BoxSize*=self.ScaleFactor
        if kwargs.get('convert_to_comoving')!=None:
            self.convert_to_comoving=kwargs.get('convert_to_comoving')
            if self.convert_to_comoving==True:
                self.pos/=self.ScaleFactor
                self.BoxSize/=self.ScaleFactor
        if kwargs.get('convert_to_per_littleh')!=None:
            self.convert_to_per_littleh=kwargs.get('convert_to_per_littleh')
            if self.convert_to_per_littleh==True:
                self.pos*=self.HubbleParam
                self.mass*=self.HubbleParam
                self.BoxSize*=self.HubbleParam
        if kwargs.get('convert_to_littleh')!=None:
            self.convert_to_littleh=kwargs.get('convert_to_littleh')
            if self.convert_to_littleh==True:
                self.pos/=self.HubbleParam
                self.mass/=self.HubbleParam
                self.BoxSize/=self.HubbleParam
                
def select_particles(val,valoffset,size,geometry,**kwargs):
    dval=val-valoffset
    # First check for periodicity
    if kwargs.get('periodic')==True and kwargs.get('scale_length')!=None:
        scale_length=kwargs.get('scale_length')
        dval=np.where(dval>0.5*scale_length,dval-scale_length,dval)
        dval=np.where(dval<-0.5*scale_length,dval+scale_length,dval)
    else:
        print('Ignoring periodicity')
    # Impose cut based on desired geometry
    if geometry=='cubic':
        ipick=np.logical_and(np.abs(dval[:,0])<size,np.abs(dval[:,1])<size)
        ipick=np.logical_and(ipick,np.abs(dval[:,2])<size)
    elif geometry=='spherical':
        r2=dval[:,0]**2+dval[:,1]**2+dval[:,2]**2
        ipick=np.where(r2<size*size)[0]
    else:
        print('Undefined geometry')
        ipick=None
    return ipick

    
