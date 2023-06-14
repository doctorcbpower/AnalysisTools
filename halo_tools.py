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
    def __init__(self,halocatfilename,halocatfileformat,comoving_units=True,**kwargs):
        self.halocatfilename=halocatfilename
        self.halocatfileformat=halocatfileformat
        self.comoving_units=comoving_units                

    def ReadHaloCatalogue(self):
        if self.halocatfileformat=='SubFind':
            with h5py.File(self.halocatfilename,'r') as f:
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
                SubhaloRankInGr=f['Subhalo/SubhaloRankInGr'][()]  # 0 is a field halo
                
                self.SubhaloGroupNr=f['Subhalo/SubhaloGroupNr'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubhaloLen=f['Subhalo/SubhaloLen'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubhaloOffsetType=f['Subhalo/SubhaloOffsetType'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubPos=f['Subhalo/SubhaloPos'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubVel=f['Subhalo/SubhaloVel'][()]#[np.where(SubhaloRankInGr!=0)[0]]
                self.SubMass=f['Subhalo/SubhaloMass'][()]#[np.where(SubhaloRankInGr!=0)[0]]

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
           if os.path.exists(filename==False):
               filename+='.0'
               
           with h5py.File(filename,'r') as f:               
                # # Read Halo Properties
               Structuretype=f['Structuretype'][()]  # 10 is a field halo
               HubbleParam=np.float32(f['SimulationInfo'].attrs['h_val'])
               ScaleFactor=np.float32(f['SimulationInfo'].attrs['ScaleFactor'])
               self.GroupAscale=ScaleFactor*np.ones(len(np.where(Structuretype==10)[0]))
               self.GroupNsubs=f['numSubStruct'][()][np.where(Structuretype==10)[0]]
               self.GroupMass=f['Mass_tot'][()][np.where(Structuretype==10)[0]]               
               self.GroupPos=(np.array([f['Xcminpot'][()],f['Ycminpot'][()],f['Zcminpot'][()]]).T)[np.where(Structuretype==10)[0]]
               self.GroupVel=(np.array([f['VXcminpot'][()],f['VYcminpot'][()],f['VZcminpot'][()]]).T)[np.where(Structuretype==10)[0]]         
               self.GroupM200=f['Mass_200crit'][()][np.where(Structuretype==10)[0]]
               self.GroupR200=f['R_200crit'][()][np.where(Structuretype==10)[0]]
               self.GroupLen=f['npart'][()][np.where(Structuretype==10)[0]]
               if self.comoving_units==True:
                   self.GroupPos*=HubbleParam
                   self.GroupPos/=ScaleFactor
                   self.GroupR200*=HubbleParam
                   self.GroupR200/=ScaleFactor
                   self.GroupM200*=HubbleParam
                   self.GroupMass*=HubbleParam
                   
               # # Read Subhalo Properties
               HostID=f['hostHaloID'][()]  # -1 is a field halo
               HostIDOffset=np.min(HostID[np.where(HostID!=-1)[0]])
               self.SubhaloGroupNr=(HostID-HostIDOffset)[np.where(Structuretype!=10)[0]]
               self.SubhaloLen=f['npart'][()][np.where(Structuretype!=10)[0]]
               self.SubPos=(np.array([f['Xcminpot'][()],f['Ycminpot'][()],f['Zcminpot'][()]]).T)[np.where(Structuretype!=10)[0]]
               self.SubVel=(np.array([f['VXcminpot'][()],f['VYcminpot'][()],f['VZcminpot'][()]]).T)[np.where(Structuretype!=10)[0]]         
               self.SubMass=f['Mass_tot'][()][np.where(Structuretype!=10)[0]]               
               self.SubhaloLen=f['npart'][()][np.where(Structuretype!=10)[0]]               
               if self.comoving_units==True:
                   self.SubPos*=HubbleParam
                   self.SubPos/=ScaleFactor
                   self.SubMass*=HubbleParam
