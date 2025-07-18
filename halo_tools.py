#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 13:55:13 2023

@author: cpower

Script to read and analyse halo catalogues generated using AHF, SubFind, VELOCIraptor.
"""
import h5py
import numpy as np
import os

class HaloTools:
    def __init__(self,halocatfilename,halocatfileformat,comoving_units=False,**kwargs):
        self.halocatfilename=halocatfilename
        self.halocatfileformat=halocatfileformat
        self.comoving_units=comoving_units
        self.usehalocatonly=False
        self.usesubstructure_file=False
        if 'usehalocatonly' in kwargs:
            self.usehalocatonly=kwargs['usehalocatonly']
        if 'usesubstructure_file' in kwargs:
            self.usesubstructure_file=kwargs['usesubstructure_file']


    def ReadHaloCatalogue(self):
        if self.halocatfileformat=='SubFind':
            with h5py.File(self.halocatfilename,'r') as f:
                NumFiles=f['Header'].attrs['NumFiles']
                self.BoxSize=f['Header'].attrs['BoxSize']
                self.HubbleParam=f['Parameters'].attrs['HubbleParam']

                self.TotNgroups=f['Header'].attrs['Ngroups_Total'][()]
<<<<<<< HEAD
#                print(type(f['Header'].attrs.keys()))
                if 'subgroups' in f['Header'].attrs.keys():
=======
                if any('subgroups' in x for x in list(f['Header'].attrs.keys()))==True:
>>>>>>> 5bb954d (Synchronising after too long - will start to tidy up scripts and document them)
                    self.TotNsubgroups=f['Header'].attrs['Nsubgroups_Total'][()]
                else:
                    self.TotNsubgroups=f['Header'].attrs['Nsubhalos_Total'][()]
                print('Reading data for %d groups and %d subgroups'%(self.TotNgroups,self.TotNsubgroups))
                if NumFiles>1:
                    print('Data is split across %d files'%NumFiles)

                self.GroupID=np.arange(self.TotNgroups,dtype=np.uint64)
                # # Read Halo Properties
                #self.GroupAscale=f['Group/GroupAscale'][()]
                self.GroupFirstSub=f['Group/GroupFirstSub'][()]
                self.GroupMass=f['Group/GroupMass'][()]  # This is the total mass in the FOF group
                self.GroupNsubs=f['Group/GroupNsubs'][()]
                self.GroupPos=f['Group/GroupPos'][()]
                self.GroupVel=f['Group/GroupVel'][()]
                self.GroupM200=f['Group/Group_M_Crit200'][()]
                self.GroupR200=f['Group/Group_R_Crit200'][()]
                self.GroupLen=f['Group/GroupLen'][()]
                if any('GroupOffsetType' in x for x in list(f['Group'].keys()))==True:
                    self.GroupOffsetType=f['Group/GroupOffsetType'][()]

                # # Read Subhalo Properties
                #self.SubhaloRankInGr=f['Subhalo/SubhaloRankInGr'][()]  # 0 is a field halo
<<<<<<< HEAD
                self.SubhaloGroupNr=f['Subhalo/SubhaloGroupNr'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubhaloLen=f['Subhalo/SubhaloLen'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubhaloOffsetType=f['Subhalo/SubhaloOffsetType'][()]#[np.where(SubhaloRankInGr!=0)[0]]
=======
                if any('SubhaloGroupNr' in x for x in list(f['Subhalo'].keys()))==True:
                    self.SubhaloGroupNr=f['Subhalo/SubhaloGroupNr'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                else:
                    self.SubhaloGroupNr=f['Subhalo/SubhaloGrNr'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubhaloLen=f['Subhalo/SubhaloLen'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                if any('SubhaloOffsetType' in x for x in list(f['Subhalo'].keys()))==True:
                    self.SubhaloOffsetType=f['Subhalo/SubhaloOffsetType'][()]#[np.where(SubhaloRankInGr!=0)[0]]
>>>>>>> 5bb954d (Synchronising after too long - will start to tidy up scripts and document them)
                self.SubhaloPos=f['Subhalo/SubhaloPos'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubhaloVel=f['Subhalo/SubhaloVel'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubhaloMass=f['Subhalo/SubhaloMass'][()]#[np.where(SubhaloRankInGr!=0)[0]]

        elif self.halocatfileformat=='SubFind-EAGLE':
            filename=self.halocatfilename+'.hdf5'
            if os.path.exists(filename)==False:
                filename=self.halocatfilename+'.0.hdf5'            
             
            with h5py.File(filename,'r') as f:
                NumFiles=f['Header'].attrs['NumFilesPerSnapshot']
                self.TotNgroups=f['Header'].attrs['TotNgroups'][()]
                self.TotNsubgroups=f['Header'].attrs['TotNsubgroups'][()]
                print('Reading data for %d groups and %d subgroups'%(self.TotNgroups,self.TotNsubgroups))
                if NumFiles>1:
                    print('Data is split across %d files'%NumFiles)

            self.GroupID=np.arange(self.TotNgroups,dtype=np.uint64)

            if NumFiles>1:                
                self.GroupFirstSub=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.GroupMass=np.ndarray(shape=(self.TotNgroups))
                self.GroupNsubs=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.GroupPos=np.ndarray(shape=(self.TotNgroups,3))
                self.GroupM200=np.ndarray(shape=(self.TotNgroups))
                self.GroupR200=np.ndarray(shape=(self.TotNgroups))
                self.GroupLen=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.GroupOffsetType=np.ndarray(shape=(self.TotNgroups),dtype=np.uint64)
                
                self.SubhaloGroupNr=np.ndarray(shape=(self.TotNsubgroups),dtype=np.uint32)
                self.SubhaloLen=np.ndarray(shape=(self.TotNsubgroups,6),dtype=np.uint32)
                self.SubhaloOffsetType=np.ndarray(shape=(self.TotNsubgroups),dtype=np.uint64)
                self.SubPos=np.ndarray(shape=(self.TotNsubgroups,3))
                self.SubVel=np.ndarray(shape=(self.TotNsubgroups,3))            
                self.SubMass=np.ndarray(shape=(self.TotNsubgroups))        
                    
                igstart=0
                isstart=0
                for i in range(NumFiles):
                    filename=self.halocatfilename+'.%d.hdf5'%i  
                    with h5py.File(filename,'r') as f:
                        Nsubgroups=f['Header'].attrs['Nsubgroups']
                        Ngroups=f['Header'].attrs['Ngroups']
                        igfinish=igstart+Ngroups
                        self.GroupFirstSub[igstart:igfinish]=f['FOF/FirstSubhaloID'][()]
                        self.GroupMass[igstart:igfinish]=f['FOF/GroupMass'][()]  # This is the total mass in the FOF group
                        self.GroupNsubs[igstart:igfinish]=f['FOF/NumOfSubhalos'][()]
                        self.GroupPos[igstart:igfinish]=f['FOF/GroupCentreOfPotential'][()]
                        # self.GroupVel=f['FOF/GroupVel'][()]
                        self.GroupM200[igstart:igfinish]=f['FOF/Group_M_Crit200'][()]
                        self.GroupR200[igstart:igfinish]=f['FOF/Group_R_Crit200'][()]
                        self.GroupLen[igstart:igfinish]=f['FOF/GroupLength'][()]
                        self.GroupOffsetType[igstart:igfinish]=f['FOF/GroupOffset'][()]
                        igstart=igfinish
                        
                        isfinish=isstart+Nsubgroups
                        self.SubhaloGroupNr[isstart:isfinish]=f['Subhalo/GroupNumber'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                        self.SubhaloLen[isstart:isfinish]=f['Subhalo/SubLengthType'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                        self.SubhaloOffsetType[isstart:isfinish]=f['Subhalo/SubOffset'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                        self.SubPos[isstart:isfinish]=f['Subhalo/CentreOfPotential'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                        self.SubVel[isstart:isfinish]=f['Subhalo/Velocity'][()]#[np.where(SubhaloRankInGr!=0)[0]]                
                        self.SubMass[isstart:isfinish]=f['Subhalo/Mass'][()]   
                        isstart=isfinish
            else:   
                    # # Read Halo Properties
                    # self.GroupAscale=f['FOF/GroupAscale'][()]
                    self.GroupFirstSub=f['FOF/FirstSubhaloID'][()]
                    self.GroupMass=f['FOF/GroupMass'][()]  # This is the total mass in the FOF group
                    self.GroupNsubs=f['FOF/NumOfSubhalos'][()]
                    self.GroupPos=f['FOF/GroupCentreOfPotential'][()]
                    # self.GroupVel=f['FOF/GroupVel'][()]
                    self.GroupM200=f['FOF/Group_M_Crit200'][()]
                    self.GroupR200=f['FOF/Group_R_Crit200'][()]
                    self.GroupLen=f['FOF/GroupLength'][()]
                    self.GroupOffsetType=f['FOF/GroupOffset'][()]

                   # # Read Subhalo Properties
                    SubhaloRankInGr=f['Subhalo/SubGroupNumber'][()]  # 0 is a field halo
               
                    self.SubhaloGroupNr=f['Subhalo/GroupNumber'][()][np.where(SubhaloRankInGr!=0)[0]]
                    self.SubhaloLen=f['Subhalo/SubLengthType'][()][np.where(SubhaloRankInGr!=0)[0]]
                    self.SubhaloOffsetType=f['Subhalo/SubOffset'][()][np.where(SubhaloRankInGr!=0)[0]]
                    self.SubPos=f['Subhalo/CentreOfPotential'][()][np.where(SubhaloRankInGr!=0)[0]]
                    self.SubVel=f['Subhalo/Velocity'][()][np.where(SubhaloRankInGr!=0)[0]]                
                    self.SubMass=f['Subhalo/Mass'][()][np.where(SubhaloRankInGr!=0)[0]]

        elif self.halocatfileformat=='VELOCIraptor':
           ## Get halo properties
           filename=self.halocatfilename+'.properties'
           if os.path.exists(filename)==False:
               filename+='.0'
               
           with h5py.File(filename,'r') as f:
               self.ScaleFactor=np.float32(f['SimulationInfo'].attrs['ScaleFactor'])
               self.BoxSize=np.float32(f['SimulationInfo'].attrs['Period'])
               self.HubbleParam=np.float32(f['SimulationInfo'].attrs['h_val'])
               self.Omega0=np.float32(f['SimulationInfo'].attrs['Omega_m'])
               self.OmegaBar=np.float32(f['SimulationInfo'].attrs['Omega_b'])
               NumFiles=f['Num_of_files'][0]
               self.TotNgroups=f['Total_num_of_groups'][0]
               print(f'Reading data for {self.TotNgroups} groups')
               if NumFiles>1:
                    print(f'Data is split across {NumFiles} files')

           if NumFiles>1:
                self.Structuretype=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.GroupAscale=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupMass=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupNsubs=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.GroupPos=np.ndarray(shape=(self.TotNgroups,3),dtype=np.float32)
                self.GroupPosCM=np.ndarray(shape=(self.TotNgroups,3),dtype=np.float32)
                self.GroupPosMBP=np.ndarray(shape=(self.TotNgroups,3),dtype=np.float32)
                self.GroupVel=np.ndarray(shape=(self.TotNgroups,3),dtype=np.float32)
                self.GroupVelCM=np.ndarray(shape=(self.TotNgroups,3),dtype=np.float32)
                self.GroupM200=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupSOM200=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupMFOF=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupR200=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupSOR200=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupEkin=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupEpot=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupLen=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.GroupID=np.ndarray(shape=(self.TotNgroups),dtype=np.uint64)
                self.SubhaloGroupNr=np.ndarray(shape=(self.TotNgroups),dtype=np.int64)
                    
                igstart=np.uint32(0)
                isstart=np.uint32(0)
                for i in range(NumFiles):
                    filename=self.halocatfilename+'.properties.%d'%i
                    with h5py.File(filename,'r') as f:
                        Ngroups=f['Num_of_groups'][0]
                        igfinish=igstart+Ngroups

                        self.Structuretype[igstart:igfinish]=f['Structuretype'][()]
                        self.GroupNsubs[igstart:igfinish]=f['numSubStruct'][()]
                        self.GroupMass[igstart:igfinish]=f['Mass_tot'][()]
                        self.GroupM200[igstart:igfinish]=f['Mass_200crit'][()]
                        self.GroupMFOF[igstart:igfinish]=f['Mass_FOF'][()]
<<<<<<< HEAD
                        self.GroupSOM200[igstart:igfinish]=f['Mass_200crit'][()]
=======
#                        self.GroupSOM200[igstart:igfinish]=f['Mass_200crit'][()]
>>>>>>> 5bb954d (Synchronising after too long - will start to tidy up scripts and document them)
                        self.GroupPos[igstart:igfinish]=(np.array([f['Xcminpot'][()],f['Ycminpot'][()],f['Zcminpot'][()]]).T)
                        self.GroupPosMBP[igstart:igfinish]=(np.array([f['Xcmbp'][()],f['Ycmbp'][()],f['Zcmbp'][()]]).T)
                        self.GroupPosCM[igstart:igfinish]=(np.array([f['Xc'][()],f['Yc'][()],f['Zc'][()]]).T)
                        self.GroupVel[igstart:igfinish]=(np.array([f['VXcminpot'][()],f['VYcminpot'][()],f['VZcminpot'][()]]).T)
                        self.GroupVelCM[igstart:igfinish]=(np.array([f['VXc'][()],f['VYc'][()],f['VZc'][()]]).T)
                        self.GroupR200[igstart:igfinish]=f['R_200crit'][()]
<<<<<<< HEAD
                        self.GroupSOR200[igstart:igfinish]=f['R_200crit'][()]
=======
#                        self.GroupSOR200[igstart:igfinish]=f['R_200crit'][()]
>>>>>>> 5bb954d (Synchronising after too long - will start to tidy up scripts and document them)
                        
                        self.GroupEkin[igstart:igfinish]=f['Ekin'][()]
                        self.GroupEpot[igstart:igfinish]=f['Epot'][()]
                        self.GroupLen[igstart:igfinish]=f['npart'][()]
                        self.GroupID[igstart:igfinish]=f['ID'][()]
                        self.SubhaloGroupNr[igstart:igfinish]=f['hostHaloID'][()]
                        igstart=igfinish

           else:
                   with h5py.File(filename,'r') as f:
                        self.Structuretype=f['Structuretype'][()]
                        self.GroupNsubs=f['numSubStruct'][()]
                        self.GroupMass=f['Mass_tot'][()]
                        self.GroupM200=f['Mass_200crit'][()]
<<<<<<< HEAD
                        self.GroupSOM200=f['SO_Mass_200.000000_rhocrit'][()]
=======
#                        self.GroupSOM200=f['SO_Mass_200.000000_rhocrit'][()]
>>>>>>> 5bb954d (Synchronising after too long - will start to tidy up scripts and document them)
                        self.GroupMFOF=f['Mass_FOF'][()]
                        self.GroupPos=(np.array([f['Xcminpot'][()],f['Ycminpot'][()],f['Zcminpot'][()]]).T)
                        self.GroupPosMBP=(np.array([f['Xcmbp'][()],f['Ycmbp'][()],f['Zcmbp'][()]]).T)
                        self.GroupPosCM=(np.array([f['Xc'][()],f['Yc'][()],f['Zc'][()]]).T)
                        self.GroupVel=(np.array([f['VXcminpot'][()],f['VYcminpot'][()],f['VZcminpot'][()]]).T)
                        self.GroupVelCM=(np.array([f['VXc'][()],f['VYc'][()],f['VZc'][()]]).T)
                        self.GroupR200=f['R_200crit'][()]
<<<<<<< HEAD
                        self.GroupSOR200=f['SO_R_200.000000_rhocrit'][()]
=======
#                        self.GroupSOR200=f['SO_R_200.000000_rhocrit'][()]
>>>>>>> 5bb954d (Synchronising after too long - will start to tidy up scripts and document them)
                        self.GroupEkin=f['Ekin'][()]
                        self.GroupEpot=f['Epot'][()]
                        self.GroupLen=f['npart'][()]
                        self.GroupID=f['ID'][()]
                        self.SubhaloGroupNr=f['hostHaloID'][()]

           if self.usehalocatonly==False:
            
                ## Get halo catalogue groups information
                filename=self.halocatfilename+'.catalog_groups'
                if os.path.exists(filename)==False:
                    filename+='.0'

                with h5py.File(filename,'r') as f:
                    NumFiles=f['Num_of_files'][0]
                    TotNumGroups=f['Total_num_of_groups'][0]
                    print('Reading data for %d groups'%(TotNumGroups))

                    if NumFiles>1:
                        print('Catalogue data is split across %d files'%NumFiles)
                        self.GroupSize=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                        self.GroupOffset=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                        self.GroupOffsetUnbound=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                        self.FileID=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                        igstart=np.uint32(0)
                        isstart=np.uint32(0)
                        for i in range(NumFiles):
                            filename=self.halocatfilename+'.catalog_groups.%d'%i
                            with h5py.File(filename,'r') as f:
                                NumGroups=f['Num_of_groups'][0]
                                igfinish=igstart+NumGroups
                                self.GroupSize[igstart:igfinish]=f['Group_Size'][()]
                                self.GroupOffset[igstart:igfinish]=f['Offset'][()]
                                self.GroupOffsetUnbound[igstart:igfinish]=f['Offset_unbound'][()]
                                self.FileID[igstart:igfinish]=f['File_id'][()]
                                igstart=igfinish
                    else:
                        with h5py.File(filename,'r') as f:
                            self.GroupSize=f['Group_Size'][()]
                            self.GroupOffset=f['Offset'][()]
                            self.GroupOffsetUnbound=f['Offset_unbound'][()]
                            self.FileID=f['File_id'][()]

                ## Get bound halo particles
                filename=self.halocatfilename+'.catalog_particles'
                if os.path.exists(filename)==False:
                    filename+='.0'
               
                with h5py.File(filename,'r') as f:
                    NumFiles=f['Num_of_files'][0]
                    TotalNumPartInGroups=f['Total_num_of_particles_in_all_groups'][0]
                    print('Reading data for %d particles in groups'%(TotalNumPartInGroups))
                    if NumFiles>1:
                        print('Particle IDs data is split across %d files'%NumFiles)
                        self.ParticleIDsInGroups=np.ndarray(shape=(TotalNumPartInGroups),dtype=np.uint32)
                        igstart=np.uint32(0)
                        isstart=np.uint32(0)
                        for i in range(NumFiles):
                            filename=self.halocatfilename+'.catalog_particles.%d'%i
                            with h5py.File(filename,'r') as f:
                                NumPartInGroups=f['Num_of_particles_in_groups'][0]
                                igfinish=igstart+NumPartInGroups
                                self.ParticleIDsInGroups[igstart:igfinish]=f['Particle_IDs'][()]
                                igstart=igfinish
                    else:
                        with h5py.File(filename,'r') as f:
                            self.ParticleIDsInGroups=f['Particle_IDs'][()]

                ## Get unbound halo particles
                filename=self.halocatfilename+'.catalog_particles.unbound'
                if os.path.exists(filename)==False:
                    filename+='.0'
               
                with h5py.File(filename,'r') as f:
                    NumFiles=f['Num_of_files'][0]
                    TotalNumPartInGroups=f['Total_num_of_particles_in_all_groups'][0]
                    print('Reading data for %d particles in groups'%(TotalNumPartInGroups))
                    if NumFiles>1:
                        print('Particle IDs data is split across %d files'%NumFiles)
                        self.UnboundParticleIDsInGroups=np.ndarray(shape=(TotalNumPartInGroups),dtype=np.uint32)
                        igstart=np.uint32(0)
                        isstart=np.uint32(0)
                        for i in range(NumFiles):
                            filename=self.halocatfilename+'.catalog_particles.unbound.%d'%i
                            with h5py.File(filename,'r') as f:
                                NumPartInGroups=f['Num_of_particles_in_groups'][0]
                                igfinish=igstart+NumPartInGroups
                                self.UnboundParticleIDsInGroups[igstart:igfinish]=f['Particle_IDs'][()]
                                igstart=igfinish
                    else:
                        with h5py.File(filename,'r') as f:
                            self.UnboundParticleIDsInGroups=f['Particle_IDs'][()]

        elif self.halocatfileformat=='AHF':
           ## Get halo properties
           filename=self.halocatfilename+'.AHF_halos'
               
           dtype=[("haloid",np.int64),
                   ("hostHaloID",np.int64),
                   ("numSubStruct",np.int64),
                   ("mhalo",np.float32),
                   ("npart",np.int64),
                   ("xcen",np.float32),
                   ("ycen",np.float32),
                   ("zcen",np.float32),
                   ("vxcen",np.float32),
                   ("vycen",np.float32),
                   ("vzcen",np.float32),
                   ("rhalo",np.float32),
                   ("com_offset",np.float32),
                   ("vmax",np.float32)]
    
           (haloid,
            hostHaloID,
            numSubStruct,
            mhalo,
            npart,
            xcen,ycen,zcen,
            vxcen,vycen,vzcen,
            rhalo,
            ds,
            vmax
            )=np.loadtxt(filename,
                            usecols=(0,1,2,3,4,5,6,7,8,9,10,11,15,16),
                            unpack=True,
                            dtype=dtype)

           self.TotNgroups=np.where(hostHaloID==0)[0].size
           print('Reading data for %d groups'%(self.TotNgroups))

           self.GroupID=haloid
           self.GroupNsubs=numSubStruct
           self.GroupParentHaloID=hostHaloID
           self.GroupMass=mhalo
           self.GroupPos=(np.array([xcen,ycen,zcen]).T)
           self.GroupVel=(np.array([vxcen,vycen,vzcen]).T)
           self.GroupM200=mhalo
           self.GroupR200=rhalo
           self.GroupLen=npart
           
    def ReadHaloParticleIDs(self,verbose=True,**kwargs):
        if self.halocatfileformat=='AHF':
            ## Get halo properties
            filename=self.halocatfilename+'.AHF_particles'
            
            dtype=[("var1",np.uint32),("var2",np.uint64)]

            offset=0
            self.TotNGroups=np.loadtxt(filename,usecols=0,skiprows=offset,dtype=np.uint32,max_rows=1)
            
            if verbose==True:
                print('Number of groups: %d'%self.TotNGroups)

            offset+=1

            (var1,var2)=np.loadtxt(filename,usecols=(0,1),unpack=True,skiprows=offset,dtype=dtype)

            offset=0

            TotalNumPartInGroups=var1.size-self.TotNGroups

            if verbose==True:
                print('Number of particles in groups: %d'%TotalNumPartInGroups)
            
            self.GroupLen=np.ndarray(shape=(self.TotNGroups),dtype=np.uint32)
            self.GroupOffset=np.ndarray(shape=(self.TotNGroups),dtype=np.uint32)
            self.ParticleIDsInGroups=np.ndarray(shape=(TotalNumPartInGroups),dtype=np.uint32)
            self.ParticleTypeInGroups=np.ndarray(shape=(TotalNumPartInGroups),dtype=np.uint32)

            offset=0
            noffset=0
<<<<<<< HEAD
            for n in range(TotNumGroups):
                self.GroupLen[n]=var1[noffset]
                noffset+=1
                self.ParticleIDsInGroups[offset:self.GroupSize[n]+offset]=var1[noffset:self.GroupSize[n]+noffset]
                self.ParticleTypeInGroups[offset:self.GroupSize[n]+offset]=var2[noffset:self.GroupSize[n]+noffset]
=======
            for n in range(self.TotNGroups):
                self.GroupLen[n]=var1[noffset]
                noffset+=1
                self.ParticleIDsInGroups[offset:self.GroupLen[n]+offset]=var1[noffset:self.GroupLen[n]+noffset]
                self.ParticleTypeInGroups[offset:self.GroupLen[n]+offset]=var2[noffset:self.GroupLen[n]+noffset]
>>>>>>> 5bb954d (Synchronising after too long - will start to tidy up scripts and document them)
                self.GroupOffset[n]=offset
                offset+=self.GroupLen[n]
                noffset+=self.GroupLen[n]
#            dtype=[]
#
#            num_rows_to_skip=0   # Need to read in number of groups
#
#            Ngroups=np.loadtxt(filename,usecols=0,skiprows=num_rows_to_skip,dtype=np.uint32,max_rows=1)
#
#            num_rows_to_skip+=1  # Skip over number of groups
#
#            if verbose==True:
#                print('Number of groups: %d'%Ngroups)
#                        
#            num_rows_to_skip+=np.sum(self.GroupLen[:haloid])+haloid
#            
#            if verbose==True:
#                print('Skipping %d rows, %d halos'%(num_rows_to_skip,haloid))
#                
#            dtype=[("num_in_group",np.uint32),("group_id",np.uint64)]
#
#            (num_in_group,group_id)=np.loadtxt(filename,usecols=(0,1),unpack=True,skiprows=num_rows_to_skip,dtype=dtype,max_rows=1)
#            
#            if verbose==True:
#                print('Group %d, Num_in_group: %d'%(group_id,num_in_group))
#                
#            num_rows_to_skip+=1
#            num_rows_to_read=self.GroupLen[haloid]
#            
#            dtype=[("particle_id",np.uint32),("particle_type",np.int32)]
#
#            (particle_id,particle_type)=np.loadtxt(filename,usecols=(0,1),unpack=True,skiprows=num_rows_to_skip,dtype=dtype,max_rows=num_rows_to_read)
            
#            return particle_id,particle_type,num_in_group,group_id

    def GetParentHaloID(self,rank_in_catalogue,verbose=True,**kwargs):
        if self.halocatfileformat=='AHF':
            return self.GroupID[np.where(self.GroupParentHaloID==0)[0]][rank_in_catalogue]
        
    def GetAssociatedSubstructure(self,haloid,verbose=True,**kwargs):
        if self.halocatfileformat=='AHF':
            if self.usesubstructure_file==True:
                ## Get halo properties
                filename=self.halocatfilename+'.AHF_substructure'

                with open(filename) as f:
                    lines=f.read().splitlines()

                nlines=len(lines)

                grp_data=np.array([],dtype=np.int32)

                for x in lines[0:nlines:2]:
                    grp_data=np.append(grp_data,np.asarray(x.split()).astype(np.int32))

                idx=np.where(grp_data[0:nlines:2]==self.GroupID[haloid])[0][0]
                        
                if idx.size==0:
                    return -1

                if verbose==True:
                    print("GroupID: %d (%d)"%(grp_data[0:nlines:2][idx],haloid))
                    print("numSubStruct: %d"%(grp_data[1:nlines:2][idx]))

                return np.sort(np.asarray(lines[1:nlines:2][idx].split()).astype(np.int64))
            else:
                ## Get halo properties...
                if len(kwargs)==0:
                    return -1

                BoxSize=kwargs['BoxSize']
                ## First find the rank of the haloid in the '.AHF_halos' file
                idx=np.where(self.GroupID==self.GroupID[haloid])[0]
                ## Now find the group centre...
                group_cen=self.GroupPos[idx]
                ## Compute the offsets of the halos and subhalos relative to haloid
                dpos=self.GroupPos-group_cen
                
                for j in range(3):
                    dpos[:,j]=np.where(dpos[:,j]>BoxSize/2,dpos[:,j]-BoxSize,dpos[:,j])
                    dpos[:,j]=np.where(dpos[:,j]<-BoxSize/2,dpos[:,j]+BoxSize,dpos[:,j])
                
                r=np.sqrt(dpos[:,0]**2+dpos[:,1]**2+dpos[:,2]**2)
                
                halos_idx=np.where(r<self.GroupR200[idx])[0]
                                
                if verbose==True:
                    print("GroupID: %d (%d)"%(self.GroupID[idx],self.GroupID[haloid]))
                    print("numSubStruct: %d"%(halos_idx.size-1))

                return np.delete(halos_idx,np.where(halos_idx==idx)[0]),self.GroupID[np.delete(halos_idx,np.where(halos_idx==idx)[0])]
        elif self.halocatfileformat=='SubFind':
            ## Get halo properties...
            if len(kwargs)==0:
                return -1

            BoxSize=kwargs['BoxSize']
            ## First find the rank of the haloid in the '.AHF_halos' file
            idx=np.where(self.GroupID==self.GroupID[haloid])[0]
            ## Now find the group centre...
            group_cen=self.GroupPos[idx]
            ## Compute the offsets of the halos and subhalos relative to haloid
            dpos=self.SubhaloPos-group_cen
                
            for j in range(3):
                dpos[:,j]=np.where(dpos[:,j]>BoxSize/2,dpos[:,j]-BoxSize,dpos[:,j])
                dpos[:,j]=np.where(dpos[:,j]<-BoxSize/2,dpos[:,j]+BoxSize,dpos[:,j])
                
            r=np.sqrt(dpos[:,0]**2+dpos[:,1]**2+dpos[:,2]**2)
                
            halos_idx=np.where(r<self.GroupR200[idx])[0]
                                
            if verbose==True:
                print("GroupID: %d (%d)"%(self.GroupID[idx],self.GroupID[haloid]))
                print("numSubStruct: %d"%(halos_idx.size-1))

                return np.delete(halos_idx,np.where(halos_idx==idx)[0]),self.GroupID[np.delete(halos_idx,np.where(halos_idx==idx)[0])]
 
