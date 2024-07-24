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
import struct

SOLAR_MASS_IN_CGS=1.989e33
YEAR_IN_CGS=60*60*24*365
KPC_IN_CGS=3.0856e21
KM_PER_SEC_IN_CGS=1.e5

class SnapshotTools:
    '''
    Read simulation snapshots - can be in HDF5 or GADGET binary format.
    
        Parameters:
            snapfilename - path of the snapshot
            snapfileformat - format of the snapshot; can be 1 or 2 (binary)
            or 3 (HDF5), or 'SNAP1/SNAP2' (binary) of 'HDF5'
            convention - 'SWIFT', 'GADGET4', 'AREPO', otherwise assumes 'GADGET2/3'
        
        If the file is in HDF5 formay, omit the .hdf5 suffix.
    '''
    def __init__(self,snapfilename,snapfileformat,gas_type=0,dm_type=1,star_type=4,bh_type=5,**kwargs):
        self.snapfilename=snapfilename
        self.snapfileformat=snapfileformat
        self.gas_type=gas_type
        self.dm_type=dm_type
        self.star_type=star_type
        self.bh_type=bh_type
        self.convention='GADGET'
        self.positions_only=False
        self.hires_only=False
        self.get_ptypes=False
        self.extra_blocks=[]
        self.positions_type='float32'
        self.pids_type=32

        # Definition of units in CGS
        self.unit_length_in_cgs=KPC_IN_CGS  # kpc
        self.unit_mass_in_cgs=1e10*SOLAR_MASS_IN_CGS # 1e10 Msol
        self.unit_velocity_in_cgs=KM_PER_SEC_IN_CGS # km/s
        self.unit_time_in_cgs=self.unit_length_in_cgs/self.unit_velocity_in_cgs
        self.unit_density_in_cgs=self.unit_mass_in_cgs/self.unit_length_in_cgs**3
        self.unit_sfr_in_cgs=SOLAR_MASS_IN_CGS/YEAR_IN_CGS
        
        if 'positions_type' in kwargs:
            self.positions_type=kwargs.get('positions_type')
        print('Assuming positions are type %s'%self.positions_type)

        if 'pids_type' in kwargs:
            self.pids_type=kwargs.get('pids_type')
        print('Assuming particle IDs are type %d bit'%self.pids_type)
        
        if self.snapfileformat=='HDF5' or self.snapfileformat=='3':
            self.convention=kwargs.get('convention')
        if kwargs.get('positions_only'):
            self.positions_only=kwargs.get('positions_only')
        if kwargs.get('hires_only'):
            self.hires_only=kwargs.get('hires_only')
        if kwargs.get('extra_blocks'):
            self.extra_blocks=kwargs.get('extra_blocks')
        if kwargs.get('get_ptypes'):
            self.get_ptypes=kwargs.get('get_ptypes')


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
                self.NumFiles=f['Header'].attrs['NumFilesPerSnapshot']
                NameOfMassBlock='Mass'
                ii=np.where(f['Header'].attrs['MassTable'][()]==0)[0]
                for indx in ii:
                    if f['Header'].attrs['NumPart_Total'][indx]>0:
                        if 'Masses' in list(f['PartType%d'%indx].keys()):
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
                elif self.convention=='GADGET4' or self.convention=='AREPO':
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
                if self.NumFiles>1:
                    print('Data is split across %d files'%self.NumFiles)
                self.ispotential=False
                if 'Potential' in f['PartType1'].keys():
                    self.ispotential=True
   
            NumPart=np.sum(self.NumPart_Total)
            print('Number of particles: %010d'%NumPart)
            print('Number of particle types: %d'%self.NumPartType)
            idx_with_mass=np.where(self.MassTable>0)[0]
            NumPart_InMassBlock=np.sum(self.NumPart_Total[idx_with_mass])
            if NumPart_InMassBlock>0:
                print('Number of particles in mass block: %010d'%NumPart_InMassBlock)
                
            if self.hires_only==True:
                NumPart=np.sum(self.NumPart_Total[:2])
                print('Number of HIRES particles: %010d'%NumPart)

            self.pos=np.ndarray(shape=(NumPart,3))
            if self.positions_only==False:
                self.vel=np.ndarray(shape=(NumPart,3))
                if self.pids_type==32:
                    self.pids=np.ndarray(shape=(NumPart),dtype=np.uint32)
                else:
                    self.pids=np.ndarray(shape=(NumPart),dtype=np.uint64)
                self.mass=np.ndarray(shape=(NumPart))
            
            if self.NumPart_Total[0]>0 and self.positions_only==False:
                self.u=np.ndarray(shape=(self.NumPart_Total[0]))
                self.rho=np.ndarray(shape=(self.NumPart_Total[0]))
                self.smoothinglength=np.ndarray(shape=(self.NumPart_Total[0]))

            if self.get_ptypes==True:
                self.ptype=np.ones(shape=(NumPart),dtype=np.int32)
                ioffset=np.zeros(shape=(self.NumPartType+1),dtype=np.uint64)
                ioffset[1:]=np.cumsum(self.NumPart_Total)
                for i in range(self.NumPartType):
                    self.ptype[ioffset[i]:ioffset[i+1]]=self.ptype[ioffset[i]:ioffset[i+1]]*i
                
            isstellarage=False
            ismetallicity=False
            issfr=False
            isstellargens=False
            isstellarinitmass=False
            ispotential=False
            
            if len(self.extra_blocks)>0:
                print('Loading extra blocks: %s'%self.extra_blocks)
                if 'AGE' in self.extra_blocks:
                    self.stellarage=np.ndarray(shape=(self.NumPart_Total[self.star_type]))
                    isstellarage=True
                if 'Z' in self.extra_blocks:
                    self.gas_metallicity=np.ndarray(shape=(self.NumPart_Total[self.gas_type]))
                    self.stellar_metallicity=np.ndarray(shape=(self.NumPart_Total[self.star_type]))
                    ismetallicity=True
                if 'SFR' in self.extra_blocks:
                    self.gas_sfr=np.ndarray(shape=(self.NumPart_Total[self.gas_type]))
                    issfr=True
                if 'STELLARGENS' in self.extra_blocks:
                    self.stellargen=np.ndarray(shape=(self.NumPart_Total[self.star_type]),dtype=np.int32)
                    isstellargens=True
                if 'INIT_MASS' in self.extra_blocks:
                    self.stellarinitmass=np.ndarray(shape=(self.NumPart_Total[self.star_type]),dtype=np.float32)
                    isstellarinitmass=True
                if 'POT' in self.extra_blocks:
                    self.potential=np.ndarray(shape=(NumPart))
                    ispotential=True
                
            istart=np.zeros(self.NumPartType,dtype=np.uint64)
            offset=0
            
            for i in range(1,self.NumPartType):
                offset+=self.NumPart_Total[i-1]
                if self.NumPart_Total[i]>0:
                    istart[i]=offset
            ifinish=np.copy(istart)
            
            jstart=0
            jfinish=np.copy(jstart)
            
            if self.NumFiles>1:
                for i in range(self.NumFiles):
                    filename=self.snapfilename+'.%d.hdf5'%i
                    print('Reading in file %s...'%filename)
                    with h5py.File(filename,'r') as f:
                        NumPart_ThisFile=f['Header'].attrs['NumPart_ThisFile'][()]
                        for itype in range(self.NumPartType):
                            if self.hires_only==True and itype>1:
                                continue
                            if NumPart_ThisFile[itype]>0:
                                ifinish[itype]=istart[itype]+NumPart_ThisFile[itype]
                                self.pos[istart[itype]:ifinish[itype]]=f['PartType%d/Coordinates'%itype][()]
                                if self.positions_only==True:
                                    continue
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
                                else:
                                    self.mass[istart[itype]:ifinish[itype]]=self.MassTable[itype]
                                if self.ispotential==True:
                                    self.potential[istart[itype]:ifinish[itype]]=f['PartType%d/Potential'%itype][()]
                                if itype==self.gas_type:
                                    if self.convention=='SWIFT':
                                        self.u[istart[itype]:ifinish[itype]]=f['PartType%d/InternalEnergies'%itype][()]
                                        self.rho[istart[itype]:ifinish[itype]]=f['PartType%d/Densities'%itype][()]
                                    else:
                                        self.u[istart[itype]:ifinish[itype]]=f['PartType%d/InternalEnergy'%itype][()]
                                        self.rho[istart[itype]:ifinish[itype]]=f['PartType%d/Density'%itype][()]
                                        self.smoothinglength[istart[itype]:ifinish[itype]]=f['PartType%d/SmoothingLength'%itype][()]
 
                                    if ismetallicity==True:
                                        if f['PartType%d/Metallicity'%itype].ndim>1:
                                            self.gas_metallicity[istart[itype]:ifinish[itype]]=f['PartType%d/Metallicity'%itype][:,0]
                                        else:
                                            self.gas_metallicity[istart[itype]:ifinish[itype]]=f['PartType%d/Metallicity'%itype][()]
                                    if issfr==True:
                                        self.gas_sfr[istart[itype]:ifinish[itype]]=f['PartType%d/StarFormationRate'%itype][()]
                                if itype==self.star_type:
                                    jfinish=jstart+NumPart_ThisFile[itype]
                                    if ismetallicity==True:
                                        if f['PartType%d/Metallicity'%itype].ndim>1:
                                            self.stellar_metallicity[jstart:jfinish]=f['PartType%d/Metallicity'%itype][:,0]
                                        else:
                                            self.stellar_metallicity[jstart:jfinish]=f['PartType%d/Metallicity'%itype][()]
                                    if isstellarage==True:
                                        self.stellarage[jstart:jfinish]=f['PartType%d/StellarFormationTime'%itype][()]
                                    if isstellargens==True:
                                        self.stellargen[jstart:jfinish]=f['PartType%d/ID_Generations'%itype][()]
                                    if isstellarinitmass==True:
                                        self.stellarinitmass[jstart:jfinish]=f['PartType%d/StellarInitMass'%itype][()]

                                    jstart=jfinish

                                istart[itype]=ifinish[itype]
            else:
                with h5py.File(filename,'r') as f:
                    for itype in range(self.NumPartType):
                        if self.hires_only==True and itype>1:
                            continue
                        if self.NumPart_Total[itype]>0:
                            ifinish[itype]=istart[itype]+self.NumPart_Total[itype]
                            self.pos[istart[itype]:ifinish[itype]]=f['PartType%d/Coordinates'%itype][()]
                            if self.positions_only==True:
                                continue

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
                            else:
                                self.mass[istart[itype]:ifinish[itype]]=self.MassTable[itype]
                            if self.ispotential==True:
                                self.potential[istart[itype]:ifinish[itype]]=f['PartType%d/Potential'%itype][()]
                                
                            if itype==self.gas_type:
                                if self.convention=='SWIFT':
                                    self.u[istart[itype]:ifinish[itype]]=f['PartType%d/InternalEnergies'%itype][()]
                                    self.rho[istart[itype]:ifinish[itype]]=f['PartType%d/Densities'%itype][()]
                                else:
                                    self.u[istart[itype]:ifinish[itype]]=f['PartType%d/InternalEnergy'%itype][()]
                                    self.rho[istart[itype]:ifinish[itype]]=f['PartType%d/Density'%itype][()]
                                    self.smoothinglength[istart[itype]:ifinish[itype]]=f['PartType%d/SmoothingLength'%itype][()]

                                if ismetallicity==True:
                                    if f['PartType%d/Metallicity'%itype].ndim>1:
                                        self.gas_metallicity[istart[itype]:ifinish[itype]]=f['PartType%d/Metallicity'%itype][:,0]
                                    else:
                                        self.gas_metallicity[istart[itype]:ifinish[itype]]=f['PartType%d/Metallicity'%itype][()]
                                if issfr==True:
                                    self.gas_sfr[istart[itype]:ifinish[itype]]=f['PartType%d/StarFormationRate'%itype][()]
                            if itype==self.star_type:
                                jfinish=jstart+self.NumPart_Total[itype]
                                if ismetallicity==True:
                                    if f['PartType%d/Metallicity'%itype].ndim>1:
                                        self.stellar_metallicity[jstart:jfinish]=f['PartType%d/Metallicity'%itype][:,0]
                                    else:
                                        self.stellar_metallicity[jstart:jfinish]=f['PartType%d/Metallicity'%itype][()]
                                if isstellarage==True:
                                    self.stellarage[jstart:jfinish]=f['PartType%d/StellarFormationTime'%itype][()]
                                if isstellargens==True:
                                    self.stellargen[jstart:jfinish]=f['PartType%d/ID_Generations'%itype][()]
                                if isstellarinitmass==True:
                                    self.stellarinitmass[jstart:jfinish]=f['PartType%d/StellarInitMass'%itype][()]
 
                                jstart=jfinish
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
                self.NumFiles=np.fromfile(f,dtype=np.int32,count=1)[0]

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

                if self.NumFiles>1:
                    print('Data is split across %d files'%self.NumFiles)
            f.close()

            idx_with_mass=np.where(self.MassTable==0)[0]    # Want to know which species are in the mass block, so
                                                            # their MassTable entries will be zero
            NumPart=np.sum(self.NumPart_Total)
            print('Number of particles: %010d'%NumPart)
            NumPart_InMassBlock=np.sum(self.NumPart_Total[idx_with_mass])
            print('Number of particles in mass block: %010d'%NumPart_InMassBlock)
            
            self.pos=np.ndarray(shape=(NumPart,3))
            self.vel=np.ndarray(shape=(NumPart,3))
            if self.pids_type==32:
                self.pids=np.ndarray(shape=(NumPart),dtype=np.uint32)
            else:
                self.pids=np.ndarray(shape=(NumPart),dtype=np.uint64)
            self.mass=np.ndarray(shape=(NumPart))
            
            if self.NumPart_Total[0]>0:
                self.u=np.ndarray(shape=(self.NumPart_Total[0]))
                self.rho=np.ndarray(shape=(self.NumPart_Total[0]))

            blocknames=self.GetBlockNames()

            isstellarage=False
            ismetallicity=False
            issfr=False
            ispotential=False
            
            if len(self.extra_blocks)>0:
                print('Loading extra blocks: %s'%self.extra_blocks)
                if 'AGE' in self.extra_blocks:
                    self.stellarage=np.ndarray(shape=(self.NumPart_Total[self.star_type]))
                    isstellarage=True
                if 'Z' in self.extra_blocks:
                    self.gas_metallicity=np.ndarray(shape=(self.NumPart_Total[self.gas_type]))
                    self.stellar_metallicity=np.ndarray(shape=(self.NumPart_Total[self.star_type]))
                    ismetallicity=True
                if 'SFR' in self.extra_blocks:
                    self.gas_sfr=np.ndarray(shape=(self.NumPart_Total[self.gas_type]))
                    issfr=True
                if 'POT' in self.extra_blocks:
                    self.potential=np.ndarray(shape=(NumPart))
                    ispotential=True

            istart=np.zeros(self.NumPartType,dtype=np.uint64)
            offset=0

            for i in range(1,self.NumPartType):
                offset+=self.NumPart_Total[i-1]
                if self.NumPart_Total[i]>0:
                    istart[i]=offset
            ifinish=np.copy(istart)

            for ifile in range(self.NumFiles):
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
                    
                    if self.NumFiles>1:
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

                    if self.pids_type==32:
                        pids_block=np.fromfile(f,dtype=np.uint32,count=NumPartInFile)
                        num_bytes=4
                    else:
                        pids_block=np.fromfile(f,dtype=np.uint64,count=NumPartInFile)
                        num_bytes=8
                        
                   # Increment beyond the IDs block
                    offset+=NumPartInFile*num_bytes
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
                        
                        if ismetallicity==True:
                            offset=blocknames['Z']+20
                            f.seek(offset,os.SEEK_SET)
                            gas_metals_block=np.fromfile(f,dtype=np.float32,count=NumPartInThisFile[self.gas_type])
                            offset+=4*NumPartInThisFile[0]
                            f.seek(offset,os.SEEK_SET)
                            stellar_metals_block=np.fromfile(f,dtype=np.float32,count=NumPartInThisFile[self.star_type])
                            
                        if isstellarage==True:
                            offset=blocknames['AGE']+20
                            f.seek(offset,os.SEEK_SET)
                            stellarage_block=np.fromfile(f,dtype=np.float32,count=NumPartInThisFile[self.star_type])

                        if issfr==True:
                            offset=blocknames['SFR']+20
                            f.seek(offset,os.SEEK_SET)
                            gas_sfr_block=np.fromfile(f,dtype=np.float32,count=NumPartInThisFile[self.gas_type])

                    if ispotential==True:
                        offset=blocknames['POT']+20
                        f.seek(offset,os.SEEK_SET)
                        potential_block=np.fromfile(f,dtype=np.float32,count=NumPartInFile)

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
                                cstart=cfinish
                            else:
                                self.mass[istart[itype]:ifinish[itype]]=self.MassTable[itype]*np.ones(NumPartInThisFile[itype])

                            if ispotential==True:
                                self.potential[istart[itype]:ifinish[itype]]=potential_block[bstart:bfinish]

                            if self.NumPart_Total[0]>0:
                                self.u[istart[itype]:ifinish[itype]]=u_block[bstart:bfinish]
                                self.rho[istart[itype]:ifinish[itype]]=rho_block[bstart:bfinish]
                                if ismetallicity==True:
                                    self.gas_metallicity[istart[itype]:ifinish[itype]]=gas_metals_block[bstart:bfinish]
                                if issfr==True:
                                    self.gas_sfr[istart[itype]:ifinish[itype]]=gas_sfr_block[bstart:bfinish]

                            if self.NumPart_Total[self.star_type]>0:
                                if ismetallicity==True:
                                    self.stellar_metallicity[istart[itype]:ifinish[itype]]=stellar_metals_block[bstart:bfinish]
                                if isstellarage==True:
                                    self.stellarage[istart[itype]:ifinish[itype]]=stellarage_block[bstart:bfinish]
                                
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
    def GetBlockNames(self):
        fileroot=self.snapfilename
        self.blocknames={}

        if self.snapfileformat=='HDF5' or self.snapfileformat=='3':
            print('HELLO')
        else:
            if self.snapfileformat!='SNAP2':
                return
                
            filename=fileroot
            
            if os.path.exists(filename)==False:
                filename=fileroot+'.0'
    
            with open(filename, mode='rb') as f:
                fileContent = f.read()
                offset=0
                filesize_in_bytes=len(fileContent)
                for i in range(20):
                    a=struct.unpack("issssii", fileContent[offset:offset+16])
                    aaa=[aa.decode('utf-8') for aa in a[1:5]]
                    tag="".join(aaa)
                    blocksize=a[5]
                    if a[0]!=a[6]:
                        print("Error")
                    self.blocknames[tag.strip()]=offset
                    offset+=2*4+a[0]+blocksize
                    if filesize_in_bytes-offset==0:
                        break
        return self.blocknames
        
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
            if ifgas-isgas>0:
                loadgas=True
            if ifstar-isstar>0:
                loadstar=True
            if ifbh-isbh>0:
                loadbh=True
        if part_type=='gas':
            if ifgas-isgas>0:
                loadgas=True
            loaddm=False
        if part_type=='star':
            if ifstar-isstar>0:
                loadstar=True
            loaddm=False
        if part_type=='bh':
            if ifbh-isbh>0:
                loadbh=True
            loaddm=False
            
        if self.ispotential==False:
            self.potential=np.zeros(shape=(np.sum(self.NumPart_Total)),dtype=np.float32)
        
        if loadgas==True:
            self.gas=self.ParticleProperties(self.pos[isgas:ifgas],
                                        self.vel[isgas:ifgas],
                                        self.pids[isgas:ifgas],
                                        self.mass[isgas:ifgas],
                                        self.potential[isgas:ifgas],
                                        internal_energy=self.u[isgas:ifgas],
                                        density=self.rho[isgas:ifgas])

        if loaddm==True:
            self.dm=self.ParticleProperties(self.pos[isdm:ifdm],
                                        self.vel[isdm:ifdm],
                                        self.pids[isdm:ifdm],
                                        self.mass[isdm:ifdm],
                                        self.potential[isdm:ifdm])
                                        

        if loadstar==True:
            self.star=self.ParticleProperties(self.pos[isstar:ifstar],
                                        self.vel[isstar:ifstar],
                                        self.pids[isstar:ifstar],
                                        self.mass[isstar:ifstar],
                                        self.potential[isstar:ifstar])
        
        if loadbh==True:
            self.bh=self.ParticleProperties(self.pos[isbh:ifbh],
                                        self.vel[isbh:ifbh],
                                        self.pids[isbh:ifbh],
                                        self.mass[isbh:ifbh],
                                        self.potential[isbh:ifbh])

        
    class ParticleProperties:
        def __init__(self,pos,vel,pids,mass,potential,**kwargs):
            self.pos=pos
            self.vel=vel
            self.pids=pids
            self.mass=mass
            self.potential=potential
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

    def WriteSnapshot(self,output_file,sim_type,npart_type,masstable_type,idx,idx_type,**kwargs):
        '''
        Write data to a single snapshot.
        '''
        filename=output_file+'.hdf5'
        print('Writing data to %s'%filename)
    
        with h5py.File(filename,'w') as f:
            header=f.create_group('Header')
            header.attrs['NumFilesPerSnapshot']=self.NumFiles
            header.attrs['NumPart_Total']=npart_type
            header.attrs['NumPart_Total_HighWord']=np.zeros(6,dtype=np.int32)
            header.attrs['Flag_Entropy_ICs']=int(1)
            header.attrs['MassTable']=self.MassTable
            
            if sim_type=='SWIFT':
                cosmo=f.create_group('Cosmology')
                header.attrs['Scale-factor']=self.ScaleFactor
                header.attrs['BoxSize']=self.BoxSize*np.ones(3)
                header.attrs['MassTable']=self.MassTable
                header.attrs['Dimension']=3
                cosmo.attrs['Omega_cdm']=self.OmegaDM
                cosmo.attrs['Omega_b']=self.OmegaBar
                cosmo.attrs['Omega_lambda']=self.OmegaLambda
                cosmo.attrs['h']=self.HubbleParam
            elif sim_type=='GADGET4':
                params=f.create_group('Parameters')
                header.attrs['Time']=self.ScaleFactor
                header.attrs['BoxSize']=self.BoxSize
                params.attrs['Omega0']=self.Omega0
                params.attrs['OmegaLambda']=self.OmegaLambda
                params.attrs['HubbleParam']=self.HubbleParam
            else:
                header.attrs['Time']=self.ScaleFactor
                header.attrs['Redshift']=1./self.ScaleFactor-1.
                header.attrs['BoxSize']=self.BoxSize
                header.attrs['Omega0']=self.Omega0
                header.attrs['OmegaLambda']=self.OmegaLambda
                header.attrs['HubbleParam']=self.HubbleParam
                header.attrs['NumPart_ThisFile']=npart_type
                header.attrs['Flag_Cooling']=int(0)
                header.attrs['Flag_StellarAge']=int(0)
                header.attrs['Flag_Sfr']=int(0)
                header.attrs['Flag_Metals']=int(0)
                header.attrs['Flag_Feedback']=int(0)                
                header.attrs['Flag_DoublePrecision']=int(0)

            if kwargs.get('selection')!=None:
                header.attrs['Halo Centre']=kwargs.get('selection')[0:3]
                header.attrs['Halo Systemic Velocity']=kwargs.get('selection')[4:6]
                header.attrs['Halo Extent']=kwargs.get('selection')[6]
            
            header.attrs['RunLabel']=sim_type

            iperiodic=1
            if kwargs.get('periodic')!=None:
                if kwargs.get('periodic')==False:
                    iperiodic=0
            header.attrs['Periodic']=iperiodic

            NameOfMassBlock='Masses'
            if kwargs.get('NameOfMassBlock')!=None:
                NameOfMassBlock=kwargs.get('NameOfMassBlock')
                
            NumPart=np.sum(npart_type)
            NumPartType=npart_type.size
            print('Number of particles: %010d'%NumPart)
            print('Number of particle types: %d'%NumPartType)
            idx_with_mass=np.where(masstable_type>0)[0]
            NumPart_InMassBlock=np.sum(npart_type[idx_with_mass])
            if NumPart_InMassBlock>0:
                print('Number of particles in mass block: %010d'%NumPart_InMassBlock)
            
            for i in range(NumPartType):
                if npart_type[i]>0:
                    group=f.create_group('PartType%d'%i)
                    
                    # Positions block
                    data_pos=group.create_dataset('Coordinates',data=self.pos[idx][idx_type[i]:idx_type[i+1]])
                    data_pos.attrs['CGSConversionFactor']=self.unit_length_in_cgs
                    data_pos.attrs['aexp-scale-exponent']=1
                    data_pos.attrs['h-scale-exponent']=-1

                    # Velocity block
                    data_vel=group.create_dataset('Velocities',data=self.vel[idx][idx_type[i]:idx_type[i+1]])
                    data_vel.attrs['CGSConversionFactor']=self.unit_velocity_in_cgs
                    data_vel.attrs['aexp-scale-exponent']=0.5
                    data_vel.attrs['h-scale-exponent']=0

                    # Particle IDs block
                    data_pids=group.create_dataset('ParticleIDs',data=self.pids[idx][idx_type[i]:idx_type[i+1]])
                    data_pids.attrs['CGSConversionFactor']=1
                    data_pids.attrs['aexp-scale-exponent']=0
                    data_pids.attrs['h-scale-exponent']=0

                    # Mass block
                    data_mass=group.create_dataset(NameOfMassBlock,data=self.mass[idx][idx_type[i]:idx_type[i+1]])
                    data_mass.attrs['CGSConversionFactor']=self.unit_mass_in_cgs
                    data_mass.attrs['aexp-scale-exponent']=0
                    data_mass.attrs['h-scale-exponent']=-1

                    if i==0:
                        # Internal energies block
                        data_u=group.create_dataset('InternalEnergy',data=self.u[0:idx_type[1]])
                        data_u.attrs['CGSConversionFactor']=self.unit_velocity_in_cgs**2
                        data_u.attrs['aexp-scale-exponent']=0
                        data_u.attrs['h-scale-exponent']=0

                        # Density block
                        data_density=group.create_dataset('Density',data=self.rho[0:idx_type[1]])
                        data_density.attrs['CGSConversionFactor']=self.unit_density_in_cgs
                        data_density.attrs['aexp-scale-exponent']=3
                        data_density.attrs['h-scale-exponent']=2

                        # Smoothing length block
                        data_smoothinglength=group.create_dataset('SmoothingLength',data=self.smoothinglength[0:idx_type[1]])
                        data_smoothinglength.attrs['CGSConversionFactor']=self.unit_length_in_cgs
                        data_smoothinglength.attrs['aexp-scale-exponent']=1
                        data_smoothinglength.attrs['h-scale-exponent']=-1

                        # Gas metals block
                        data_gas_metals=group.create_dataset('Metallicity',data=self.gas_metallicity[0:idx_type[1]])
                        data_gas_metals.attrs['CGSConversionFactor']=1
                        data_gas_metals.attrs['aexp-scale-exponent']=0
                        data_gas_metals.attrs['h-scale-exponent']=-0

                        
                        # Star formation rate block
                        data_gas_sfr=group.create_dataset('StarFormationRate',data=self.gas_sfr[0:idx_type[1]])
                        data_gas_sfr.attrs['CGSConversionFactor']=self.unit_sfr_in_cgs
                        data_gas_sfr.attrs['aexp-scale-exponent']=0
                        data_gas_sfr.attrs['h-scale-exponent']=0

                    if i==4:
                        # Stellar metals block
                        data_stellar_metals=group.create_dataset('Metallicity',data=self.stellar_metallicity[idx_type[i]:idx_type[i+1]])
                        data_stellar_metals.attrs['CGSConversionFactor']=1
                        data_stellar_metals.attrs['aexp-scale-exponent']=0
                        data_stellar_metals.attrs['h-scale-exponent']=0

                        # Stellar age block
                        data_stellar_age=group.create_dataset('StellarFormationTime',data=self.stellarage[idx_type[i]:idx_type[i+1]])
                        data_stellar_age.attrs['CGSConversionFactor']=1
                        data_stellar_age.attrs['aexp-scale-exponent']=0
                        data_stellar_age.attrs['h-scale-exponent']=0

                        # Stellar age block
                        data_stellar_initmass=group.create_dataset('StellarInitMass',data=self.stellarinitmass[idx_type[i]:idx_type[i+1]])
                        data_stellar_initmass.attrs['CGSConversionFactor']=self.unit_mass_in_cgs
                        data_stellar_initmass.attrs['aexp-scale-exponent']=0
                        data_stellar_initmass.attrs['h-scale-exponent']=-1


## Assume that size is the boxsize if cubic, radius if spherical.
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
    if kwargs.get('return_ptype')==True:
        return ipick,kwargs.get('ptype')[ipick]
    else:
        return ipick
    
def place_points_in_mesh(pos,pos_offset,size,mesh_dimension,**kwargs):
    return np.fix(mesh_dimension*(pos-pos_offset)/size)
    
