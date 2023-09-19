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

class HaloTools:
    def __init__(self,halocatfilename,halocatfileformat,comoving_units=False,**kwargs):
        self.halocatfilename=halocatfilename
        self.halocatfileformat=halocatfileformat
        self.comoving_units=comoving_units                

    def ReadHaloCatalogue(self):
        if self.halocatfileformat=='SubFind':
            with h5py.File(self.halocatfilename,'r') as f:
                print(f.keys())
                self.BoxSize=f['Header'].attrs['BoxSize']
                self.HubbleParam=f['Parameters'].attrs['HubbleParam']

                # # Read Halo Properties
                self.GroupAscale=f['Group/GroupAscale'][()]
                self.GroupFirstSub=f['Group/GroupFirstSub'][()]
                self.GroupMass=f['Group/GroupMass'][()]  # This is the total mass in the FOF group
                self.GroupNsubs=f['Group/GroupNsubs'][()]
                self.GroupPos=f['Group/GroupPos'][()]
                self.GroupVel=f['Group/GroupVel'][()]
                self.GroupM200=f['Group/Group_M_Crit200'][()]
                self.GroupR200=f['Group/Group_R_Crit200'][()]
                self.GroupLen=f['Group/GroupLen'][()]
                self.GroupOffsetType=f['Group/GroupOffsetType'][()]

                # # Read Subhalo Properties
                self.SubhaloRankInGr=f['Subhalo/SubhaloRankInGr'][()]  # 0 is a field halo
                self.SubhaloGroupNr=f['Subhalo/SubhaloGroupNr'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubhaloLen=f['Subhalo/SubhaloLen'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubhaloOffsetType=f['Subhalo/SubhaloOffsetType'][()]#[np.where(SubhaloRankInGr!=0)[0]]
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
#               print(list(f.keys()))
               self.ScaleFactor=np.float32(f['SimulationInfo'].attrs['ScaleFactor'])
               self.BoxSize=np.float32(f['SimulationInfo'].attrs['Period'])
               self.HubbleParam=np.float32(f['SimulationInfo'].attrs['h_val'])
               self.Omega0=np.float32(f['SimulationInfo'].attrs['Omega_m'])
               self.OmegaBar=np.float32(f['SimulationInfo'].attrs['Omega_b'])
               NumFiles=f['Num_of_files'][0]
               self.TotNgroups=f['Total_num_of_groups'][0]
               print('Reading data for %d groups'%(self.TotNgroups))
               if NumFiles>1:
                    print('Data is split across %d files'%NumFiles)

           if NumFiles>1:
#                self.GroupFirstSub=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.Structuretype=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.GroupAscale=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupMass=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupNsubs=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.GroupPos=np.ndarray(shape=(self.TotNgroups,3),dtype=np.float32)
                self.GroupVel=np.ndarray(shape=(self.TotNgroups,3),dtype=np.float32)
                self.GroupM200=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupR200=np.ndarray(shape=(self.TotNgroups),dtype=np.float32)
                self.GroupLen=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
                self.GroupID=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)

#                self.GroupOffsetType=np.ndarray(shape=(self.TotNgroups),dtype=np.uint64)
                
                self.SubhaloGroupNr=np.ndarray(shape=(self.TotNgroups),dtype=np.uint32)
#                self.SubhaloLen=np.ndarray(shape=(self.TotNsubgroups,6),dtype=np.uint32)
#                self.SubhaloOffsetType=np.ndarray(shape=(self.TotNsubgroups),dtype=np.uint64)
#                self.SubPos=np.ndarray(shape=(self.TotNsubgroups,3))
#                self.SubVel=np.ndarray(shape=(self.TotNsubgroups,3))
#                self.SubMass=np.ndarray(shape=(self.TotNsubgroups))
                    
                igstart=np.uint32(0)
                isstart=np.uint32(0)
                for i in range(NumFiles):
                    filename=self.halocatfilename+'.properties.%d'%i
                    with h5py.File(filename,'r') as f:
#                        print(list(f.keys()))
#                        Nsubgroups=f['Header'].attrs['Nsubgroups']
                        Ngroups=f['Num_of_groups'][0]
                        igfinish=igstart+Ngroups

#                        self.GroupFirstSub[igstart:igfinish]=f['FOF/FirstSubhaloID'][()]
#                        self.GroupAscale[igstart:igfinish]
                        self.Structuretype[igstart:igfinish]=f['Structuretype'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupNsubs[igstart:igfinish]=f['numSubStruct'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupMass[igstart:igfinish]=f['Mass_tot'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupPos[igstart:igfinish]=(np.array([f['Xcminpot'][()],f['Ycminpot'][()],f['Zcminpot'][()]]).T)#[np.where(Structuretype==10)[0]]
                        self.GroupVel[igstart:igfinish]=(np.array([f['VXcminpot'][()],f['VYcminpot'][()],f['VZcminpot'][()]]).T)#[np.where(Structuretype==10)[0]]
                        self.GroupM200[igstart:igfinish]=f['Mass_200crit'][()]#[np.where(Structuretype==10)[0]]
#                        self.GroupMFOF[igstart:igfinish]=f['Mass_FOF'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupR200[igstart:igfinish]=f['R_200crit'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupLen[igstart:igfinish]=f['npart'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupID[igstart:igfinish]=f['ID'][()]#[np.where(Structuretype==10)[0]]

#                        # self.GroupVel=f['FOF/GroupVel'][()]
#                        self.GroupM200[igstart:igfinish]=f['FOF/Group_M_Crit200'][()]
#                        self.GroupR200[igstart:igfinish]=f['FOF/Group_R_Crit200'][()]
#                        self.GroupLen[igstart:igfinish]=f['FOF/GroupLength'][()]
#                        self.GroupOffsetType[igstart:igfinish]=f['FOF/GroupOffset'][()]
#                        igstart=igfinish
#
#                        isfinish=isstart+Nsubgroups
                        self.SubhaloGroupNr[igstart:igfinish]=f['hostHaloID'][()]#[np.where(SubhaloRankInGr!=0)[0]]
#                        self.SubhaloLen[isstart:isfinish]=f['Subhalo/SubLengthType'][()]#[np.where(SubhaloRankInGr!=0)[0]]
#                        self.SubhaloOffsetType[isstart:isfinish]=f['Subhalo/SubOffset'][()]#[np.where(SubhaloRankInGr!=0)[0]]
#                        self.SubPos[isstart:isfinish]=f['Subhalo/CentreOfPotential'][()]#[np.where(SubhaloRankInGr!=0)[0]]
#                        self.SubVel[isstart:isfinish]=f['Subhalo/Velocity'][()]#[np.where(SubhaloRankInGr!=0)[0]]
#                        self.SubMass[isstart:isfinish]=f['Subhalo/Mass'][()]
#                        isstart=isfinish
                        igstart=igfinish

           else:
                   with h5py.File(filename,'r') as f:
#                        print(list(f.keys()))
#                        Nsubgroups=f['Header'].attrs['Nsubgroups']
                        self.Structuretype=f['Structuretype'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupNsubs=f['numSubStruct'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupMass=f['Mass_tot'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupPos=(np.array([f['Xcminpot'][()],f['Ycminpot'][()],f['Zcminpot'][()]]).T)#[np.where(Structuretype==10)[0]]
                        self.GroupVel=(np.array([f['VXcminpot'][()],f['VYcminpot'][()],f['VZcminpot'][()]]).T)#[np.where(Structuretype==10)[0]]
                        self.GroupM200=f['Mass_200crit'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupMFOF=f['Mass_FOF'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupR200=f['R_200crit'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupLen=f['npart'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupID=f['ID'][()]#[np.where(Structuretype==10)[0]]
                        self.SubhaloGroupNr=f['hostHaloID'][()]#[np.where(SubhaloRankInGr!=0)[0]]

            
#           with h5py.File(filename,'r') as f:
#                # # Read Halo Properties
#               Structuretype=f['Structuretype'][()]  # 10 is a field halo
#               HubbleParam=np.float32(f['SimulationInfo'].attrs['h_val'])
#               ScaleFactor=np.float32(f['SimulationInfo'].attrs['ScaleFactor'])
#               self.GroupAscale=ScaleFactor*np.ones(len(np.where(Structuretype==10)[0]))
#               self.GroupNsubs=f['numSubStruct'][()][np.where(Structuretype==10)[0]]
#               self.GroupMass=f['Mass_tot'][()][np.where(Structuretype==10)[0]]
#               self.GroupPos=(np.array([f['Xcminpot'][()],f['Ycminpot'][()],f['Zcminpot'][()]]).T)[np.where(Structuretype==10)[0]]
#               self.GroupVel=(np.array([f['VXcminpot'][()],f['VYcminpot'][()],f['VZcminpot'][()]]).T)[np.where(Structuretype==10)[0]]
#               self.GroupM200=f['Mass_200crit'][()][np.where(Structuretype==10)[0]]
#               self.GroupR200=f['R_200crit'][()][np.where(Structuretype==10)[0]]
#               self.GroupLen=f['npart'][()][np.where(Structuretype==10)[0]]
#               if self.comoving_units==True:
#                   self.GroupPos*=HubbleParam
#                   self.GroupPos/=ScaleFactor
#                   self.GroupR200*=HubbleParam
#                   self.GroupR200/=ScaleFactor
#                   self.GroupM200*=HubbleParam
#                   self.GroupMass*=HubbleParam
#
#               # # Read Subhalo Properties
#               HostID=f['hostHaloID'][()]  # -1 is a field halo
#               HostIDOffset=np.min(HostID[np.where(HostID!=-1)[0]])
#               self.SubhaloGroupNr=(HostID-HostIDOffset)[np.where(Structuretype!=10)[0]]
#               self.SubhaloLen=f['npart'][()][np.where(Structuretype!=10)[0]]
#               self.SubPos=(np.array([f['Xcminpot'][()],f['Ycminpot'][()],f['Zcminpot'][()]]).T)[np.where(Structuretype!=10)[0]]
#               self.SubVel=(np.array([f['VXcminpot'][()],f['VYcminpot'][()],f['VZcminpot'][()]]).T)[np.where(Structuretype!=10)[0]]
#               self.SubMass=f['Mass_tot'][()][np.where(Structuretype!=10)[0]]
#               self.SubhaloLen=f['npart'][()][np.where(Structuretype!=10)[0]]
#               if self.comoving_units==True:
#                   self.SubPos*=HubbleParam
#                   self.SubPos/=ScaleFactor
#                   self.SubMass*=HubbleParam

           ## Get halo catalogue information

           filename=self.halocatfilename+'.catalog_groups'
           if os.path.exists(filename)==False:
               filename+='.0'
               
           with h5py.File(filename,'r') as f:
               print(list(f.keys()))
               NumFiles=f['Num_of_files'][0]
               TotNumGroups=f['Total_num_of_groups'][0]
               print('Reading data for %d groups'%(TotNumGroups))
               if NumFiles>1:
                    print('Catalogue data is split across %d files'%NumFiles)

           if NumFiles>1:
                self.NumSubInGroup=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                self.GroupSize=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                self.GroupOffset=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                self.GroupOffsetUnbound=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                self.ParentHaloID=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                self.FileID=np.ndarray(shape=(TotNumGroups),dtype=np.uint32)
                igstart=np.uint32(0)
                isstart=np.uint32(0)
                for i in range(NumFiles):
                    filename=self.halocatfilename+'.catalog_groups.%d'%i
                    with h5py.File(filename,'r') as f:
                        NumGroups=f['Num_of_groups'][0]
                        igfinish=igstart+NumGroups
                        self.NumSubInGroup[igstart:igfinish]=f['Number_of_substructures_in_halo'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupSize[igstart:igfinish]=f['Group_Size'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupOffset[igstart:igfinish]=f['Offset'][()]#[np.where(Structuretype==10)[0]]
                        self.GroupOffsetUnbound[igstart:igfinish]=f['Offset_unbound'][()]#[np.where(Structuretype==10)[0]]
                        self.ParentHaloID[igstart:igfinish]=f['Parent_halo_ID'][()]#[np.where(Structuretype==10)[0]]
                        self.FileID[igstart:igfinish]=f['File_id'][()]#[np.where(Structuretype==10)[0]]
                        igstart=igfinish
           else:
                with h5py.File(filename,'r') as f:
                    self.NumSubInGroup=f['Number_of_substructures_in_halo'][()]#[np.where(Structuretype==10)[0]]
                    self.GroupSize=f['Group_Size'][()]#[np.where(Structuretype==10)[0]]
                    self.GroupOffset=f['Offset'][()]#[np.where(Structuretype==10)[0]]
                    self.GroupOffsetUnbound=f['Offset_unbound'][()]#[np.where(Structuretype==10)[0]]
                    self.ParentHaloID=f['Parent_halo_ID'][()]#[np.where(Structuretype==10)[0]]
                    self.FileID=f['File_id'][()]#[np.where(Structuretype==10)[0]]

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

           if NumFiles>1:
                self.ParticleIDsInGroups=np.ndarray(shape=(TotalNumPartInGroups),dtype=np.uint32)
                igstart=np.uint32(0)
                isstart=np.uint32(0)
                for i in range(NumFiles):
                    filename=self.halocatfilename+'.catalog_particles.%d'%i
                    with h5py.File(filename,'r') as f:
                        NumPartInGroups=f['Num_of_particles_in_groups'][0]
                        igfinish=igstart+NumPartInGroups
                        self.ParticleIDsInGroups[igstart:igfinish]=f['Particle_IDs'][()]#[np.where(Structuretype==10)[0]]
                        igstart=igfinish
           else:
                with h5py.File(filename,'r') as f:
                    self.ParticleIDsInGroups=f['Particle_IDs'][()]#[np.where(Structuretype==10)[0]]

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

           if NumFiles>1:
                self.UnboundParticleIDsInGroups=np.ndarray(shape=(TotalNumPartInGroups),dtype=np.uint32)
                igstart=np.uint32(0)
                isstart=np.uint32(0)
                for i in range(NumFiles):
                    filename=self.halocatfilename+'.catalog_particles.unbound.%d'%i
                    with h5py.File(filename,'r') as f:
                        NumPartInGroups=f['Num_of_particles_in_groups'][0]
                        igfinish=igstart+NumPartInGroups
                        self.UnboundParticleIDsInGroups[igstart:igfinish]=f['Particle_IDs'][()]#[np.where(Structuretype==10)[0]]
                        igstart=igfinish
           else:
                with h5py.File(filename,'r') as f:
                    self.UnboundParticleIDsInGroups=f['Particle_IDs'][()]#[np.where(Structuretype==10)[0]]


        elif self.halocatfileformat=='AHF':
           ## Get halo properties
           filename=self.halocatfilename+'.AHF_halos'
           print(filename)
               
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
           haloid=np.int64(haloid-1e13)
                # # Read Halo Properties
#               Structuretype=f['Structuretype'][()]  # 10 is a field halo
#               HubbleParam=np.float32(f['SimulationInfo'].attrs['h_val'])
#               ScaleFactor=np.float32(f['SimulationInfo'].attrs['ScaleFactor'])
#               self.GroupAscale=ScaleFactor*np.ones(len(np.where(Structuretype==10)[0]))
           self.GroupNsubs=numSubStruct
           self.GroupMass=mhalo
           self.GroupPos=(np.array([xcen,ycen,zcen]).T)
           self.GroupVel=(np.array([vxcen,vycen,vzcen]).T)
           self.GroupM200=mhalo
           self.GroupR200=rhalo
           self.GroupLen=npart
#               if self.comoving_units==True:
#                   self.GroupPos*=HubbleParam
#                   self.GroupPos/=ScaleFactor
#                   self.GroupR200*=HubbleParam
#                   self.GroupR200/=ScaleFactor
#                   self.GroupM200*=HubbleParam
#                   self.GroupMass*=HubbleParam
                   
#               # # Read Subhalo Properties
#               HostID=f['hostHaloID'][()]  # -1 is a field halo
#               HostIDOffset=np.min(HostID[np.where(HostID!=-1)[0]])
#               self.SubhaloGroupNr=(HostID-HostIDOffset)[np.where(Structuretype!=10)[0]]
#               self.SubhaloLen=f['npart'][()][np.where(Structuretype!=10)[0]]
#               self.SubPos=(np.array([f['Xcminpot'][()],f['Ycminpot'][()],f['Zcminpot'][()]]).T)[np.where(Structuretype!=10)[0]]
#               self.SubVel=(np.array([f['VXcminpot'][()],f['VYcminpot'][()],f['VZcminpot'][()]]).T)[np.where(Structuretype!=10)[0]]
#               self.SubMass=f['Mass_tot'][()][np.where(Structuretype!=10)[0]]
#               self.SubhaloLen=f['npart'][()][np.where(Structuretype!=10)[0]]
#               if self.comoving_units==True:
#                   self.SubPos*=HubbleParam
#                   self.SubPos/=ScaleFactor
#                   self.SubMass*=HubbleParam

    def ReadHaloParticleIDs(self,verbose=True,**kwargs):
        if self.halocatfileformat=='AHF':
            ## Get halo properties
            filename=self.halocatfilename+'.AHF_particles'
               
            dtype=[("particle_id",np.uint64),("particle_type",np.uint64)]

            num_rows_to_skip=0

            Ngroups=np.loadtxt(filename,usecols=0,skiprows=num_rows_to_skip,dtype=np.uint64,max_rows=1)

            num_rows_to_skip+=1

            if verbose==True:
                print('Number of groups: %d'%Ngroups)
                        
            if 'haloid' in kwargs:
                haloid=kwargs.get('haloid')
                num_rows_to_skip+=np.sum(self.GroupLen[:haloid-1])+haloid-1
                (num_in_group,group_id)=np.loadtxt(filename,usecols=(0,1),unpack=True,skiprows=num_rows_to_skip,dtype=np.uint64,max_rows=1)
                print(num_in_group,self.GroupLen[haloid])
                num_rows_to_skip+=1
                num_rows_to_read=self.GroupLen[haloid]
                (self.particle_id,self.particle_type)=np.loadtxt(filename,usecols=(0,1),unpack=True,dtype=dtype,skiprows=num_rows_to_skip,max_rows=num_rows_to_read)
            else:
                # Read in both columns and then slice into relevant groups based on their lengths.
                (col0,col1)=np.loadtxt(filename,usecols=(0,1),unpack=True,dtype=dtype,skiprows=1)
                print(col1[np.sum(self.GroupLen[:2-1]):np.sum(self.GroupLen[:2])+2][-1])
                
        
    def GetAssociatedSubstructure(self,haloid,verbose=True,**kwargs):
        if self.halocatfileformat=='AHF':
            ## Get halo properties
            filename=self.halocatfilename+'.AHF_substructure'

            num_rows_to_skip=0

            (group_id,num_substructure_in_group)=np.loadtxt(filename,usecols=(0,1),skiprows=num_rows_to_skip,dtype=np.uint64,max_rows=1)

            num_rows_to_skip+=1
            
            line=np.genfromtxt(filename,skiprows=num_rows_to_skip,dtype=np.uint64,max_rows=1,delimiter=' ')
            
            return line
            
    
