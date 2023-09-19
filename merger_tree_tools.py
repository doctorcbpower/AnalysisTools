#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 13:55:13 2023

@author: cpower

Script to ingest and analyse halo merger trees generated with SubFind-HBT MergerTree
and TreeFrog

"""
import h5py
import numpy as np
import os

class TreeTools:
    def __init__(self,treefilename,treefileformat,comoving_units=True,**kwargs):
        self.treefilename=treefilename
        self.treefileformat=treefileformat
        self.comoving_units=comoving_units                

    def ReadMergerTreeCatalogue(self):
        if self.treefileformat=='SubFind':
            # The SubFind-dervied tree contains a
            #   * TreeTable - with Tree IDs, Lengths, Offsets
            #   * TreeTimes - with expansion factors and redshifts
            #   * TreeHalos - with information about the subhalos, that are
            #                 the structures that are linked, and the links
            #                 between progenitors and descendents, as well
            #                 as subhalos in the same FOF structure

            with h5py.File(self.treefilename,'r') as f:

                header_info=dict(f['Header'].attrs.items())
                
                if 'Ntrees_Total' not in header_info:
                    print('Warning: no trees in file...')
                    return -1
                else:
                    print('Found %d merger trees in file...'%header_info['Ntrees_Total'])
                # Read information from the header
                self.LastSnapShotNr=header_info['LastSnapShotNr']   # Final snapshot of simulation to create tree
                self.Nhalos_ThisFile=header_info['Nhalos_ThisFile'] # Number of halos in file
                self.Nhalos_Total=header_info['Nhalos_Total']       # Number of halos in total
                self.Ntrees_ThisFile=header_info['Ntrees_ThisFile'] # Number of trees in file
                self.Ntrees_Total=header_info['Ntrees_Total']       # Number of trees in total
                self.NumFiles=header_info['NumFiles']               # Number of files

                # Read Redshift and Time
                self.Redshift=f['TreeTimes/Redshift'][()]
                self.Time=f['TreeTimes/Time'][()]

                # Read Cosmological Parameters
                self.Omega0=f['Parameters'].attrs['Omega0']
                self.OmegaLambda=f['Parameters'].attrs['OmegaLambda']
                self.OmegaBaryon=f['Parameters'].attrs['OmegaBaryon']
                self.HubbleParam=f['Parameters'].attrs['HubbleParam']
                self.BoxSize=f['Parameters'].attrs['BoxSize']
            
                # Read Halo Properties
                self.GrpNr=f['TreeHalos/GroupNr'][()]               # Group number in halo catalogue at SnapNum
                self.SubhaloNr=f['TreeHalos/SubhaloNr'][()]         # Corresponding subhalo number in halo catalogue at SnapNum
                self.SnapNum=f['TreeHalos/SnapNum'][()]             # Snapshot and halo catalogue number
                self.GrpM200=f['TreeHalos/Group_M_Crit200'][()]
                self.SubhaloMass=f['TreeHalos/SubhaloMass'][()]
                self.SubhaloPos=f['TreeHalos/SubhaloPos'][()]
                self.SubhaloVel=f['TreeHalos/SubhaloVel'][()]            
                self.SubhaloVelDisp=f['TreeHalos/SubhaloVelDisp'][()]
                self.SubhaloVmax=f['TreeHalos/SubhaloVmax'][()]
                self.SubhaloVmaxRad=f['TreeHalos/SubhaloVmaxRad'][()]
            
                # Read Tree Links
                # Descendants
                self.TreeFirstDescendant=f['TreeHalos/TreeFirstDescendant'][()]
                self.TreeDescendant=f['TreeHalos/TreeDescendant'][()]
                self.TreeNextDescendant=f['TreeHalos/TreeNextDescendant'][()]   # ID of subhalo that shares common progenitor 
                # Progenitors
                self.TreeFirstProgenitor=f['TreeHalos/TreeFirstProgenitor'][()]
                self.TreeProgenitor=f['TreeHalos/TreeProgenitor'][()]
                self.TreeNextProgenitor=f['TreeHalos/TreeNextProgenitor'][()]   # ID of subhalo that shares common descendant
                self.TreeMainProgenitor=f['TreeHalos/TreeMainProgenitor'][()]
                # FOF group components
                self.TreeFirstHaloInFOFgroup=f['TreeHalos/TreeFirstHaloInFOFgroup'][()]
                self.TreeNextHaloInFOFgroup=f['TreeHalos/TreeNextHaloInFOFgroup'][()]
                # IDs
                self.TreeHalosID=f['TreeHalos/TreeID'][()]
                self.TreeHalosIndex=f['TreeHalos/TreeIndex'][()]
                # Table Information
                # IDs
                self.TreeID=f['TreeTable/TreeID'][()]        # Unique ID of tree
                self.TreeLength=f['TreeTable/Length'][()]    # Length of tree
                self.TreeOffset=f['TreeTable/StartOffset'][()]  # Offset within catalogue of tree
        elif self.treefileformat=='TreeFrog':
            with h5py.File(self.treefilename,'r') as f:
                
                ID_Offset=f['Header/TreeBuilder'].attrs['Temporal_halo_id_value']
                HubbleParam=f['Header/Simulation'].attrs['h_val']

                self.TreeProgenitor={}
                self.TreeProgenitorSnap={}
                self.TreeHalosID={}
                self.Time={}
                self.SnapNum={}
                self.GrpM200={}
                self.SubhaloMass={}
                self.SubhaloPos={}
                self.SubhaloVel={}
                self.SubhaloVelDisp={}
                self.SubhaloVmax={}
                self.SubhaloVmaxRad={}
                
                
                for group in f.keys():                            
                    if 'Snap' in group and f[group].attrs['NHalos']>0:
                        ScaleFactor=f[group].attrs['scalefactor']
                        self.Time[group]=ScaleFactor
                        Current_ID_Offset=f[group].attrs['Snapnum']*ID_Offset+1
                        self.TreeProgenitor[group]=f[group]['Progenitor'][()]                        
                        self.TreeProgenitorSnap[group]=f[group]['ProgenitorSnap'][()]
                        self.TreeHalosID[group]=f[group]['ID'][()]
                        
                        self.SnapNum[group]=f[group].attrs['Snapnum']*np.ones(f[group].attrs['NHalos'],dtype=np.int32)
                        self.GrpM200[group]=f[group]['Mass_200crit'][()]
                        self.SubhaloMass[group]=f[group]['Mass_tot'][()]
                        self.SubhaloPos[group]=np.array([f[group]['Xcminpot'][()],f[group]['Ycminpot'][()],f[group]['Zcminpot'][()]]).T                                     
                        self.SubhaloVel[group]=np.array([f[group]['VXc'][()],f[group]['VYc'][()],f[group]['VZc'][()]]).T             
                        self.SubhaloVelDisp[group]=f[group]['sigV'][()]
                        self.SubhaloVmax[group]=f[group]['Vmax'][()]
                        self.SubhaloVmaxRad[group]=f[group]['Rmax'][()]

                        if self.comoving_units==True:
                            self.SubhaloPos[group]*=HubbleParam
                            self.SubhaloPos[group]/=ScaleFactor
                            self.SubhaloMass[group]*=HubbleParam
                         
    def TrackMainHaloProgenitor(self,MainSubhaloID,SnapshotNr):
        """ Assumes the SubFind MergerTree structure - this tracks the first progenitor of a given subhalo identified in a halo catalogue at a given snapshot number.
    
            Receives MainSubhaloID, which can be obtained from GroupFirstSub in the subhalo catalogue for a given group number; and SnapshotNr - the snapshot number of the halo catalogue in which the halo is identified.
            
            Returns Redshift, Subhalo Mass, Parent Group Mass (200 times critical), Parent ID at SnapNum, Subhalo ID at SnapNum
        """
        # Find members of SubhaloNr whose SnapNum is SnapshotNr, and the one that corresponds
        # to MainSubhaloID. Use this to identify the TreeID.
        itree=self.TreeHalosID[np.where(self.SnapNum==SnapshotNr)[0]][np.where(self.SubhaloNr[np.where(self.SnapNum==SnapshotNr)[0]]==MainSubhaloID)[0]][0]
        # Find members of SubhaloNr whose SnapNum is SnapshotNr, and the one that corresponds
        # to MainSubhaloID. Use this to identify the TreeIndex within the TreeID.
        index=self.TreeHalosIndex[np.where(self.SnapNum==SnapshotNr)[0]][np.where(self.SubhaloNr[np.where(self.SnapNum==SnapshotNr)[0]]==MainSubhaloID)[0]][0]
        # Given the TreeID, we can identify within the catalogue the offsets of the beginning and
        # end of the tree associated with TreeID
        istart=self.TreeOffset[itree]
        ifinish=self.TreeOffset[itree]+self.TreeLength[itree]
        
        # Now extract useful properties of the group and subhalo structures within tree TreeID
        grp_num=self.GrpNr[istart:ifinish]       # Group Number in catalogue at Snapshort Number
        grp_mass=self.GrpM200[istart:ifinish]    # Group M200 in catalogue at Snapshot Number
        sub_num=self.SubhaloNr[istart:ifinish]   # Subhalo Number in catalogue at Snapshot Number
        sub_mass=self.SubhaloMass[istart:ifinish]  # Subhalo Mass in catalogue at Snapshot Number
        snap_num=self.SnapNum[istart:ifinish]    # Snapshot Number
        first_prog=self.TreeMainProgenitor[istart:ifinish]  # First progenitor
        # Create arrays to contain the variables to be returned
        redshift=np.array([])
        mass=np.array([])
        m200=np.array([])
        subhalo_number=np.array([],dtype=np.uint64)
        group_number=np.array([],dtype=np.uint64)

        # Initial index is already determined at SnapNum SnapShotNr
        idx=index
        tree_length=0
        
        while first_prog[idx]!=-1:     # Loop over list until we hit the root
            redshift=np.append(redshift,self.Redshift[snap_num[idx]]) # Redshift
            mass=np.append(mass,sub_mass[idx])                        # Mass of main subhalo
            m200=np.append(m200,grp_mass[idx])                        # Mass of parent group, M200
            group_number=np.append(group_number,grp_num[idx])         # Unique ID of group at SnapNum
            subhalo_number=np.append(subhalo_number,sub_num[idx])     # Unique ID of subhalo at SnapNum
            tree_length+=1
            idx=first_prog[idx]

        zform=0.0
        
        if tree_length>10:
            zform=np.interp(0.5,np.flip(mass/mass[0]),np.flip(redshift))
        
        return redshift,mass,m200,group_number,subhalo_number,zform
            
    
    def TrackHaloDescendant(self,HaloID):
        redshift=np.array([])
        mass=np.array([])        
        index=HaloID
        while self.TreeFirstProgenitor[index] != -1:
            redshift=np.append(redshift,self.Redshift[self.SnapNum[index]])
            mass=np.append(mass,self.SubhaloMass[self.TreeFirstProgenitor[index]])
            print(index,self.GroupNr[index],self.TreeFirstDescendant[index],self.TreeNextProgenitor[index])
            index=self.TreeFirstProgenitor[index]
        return redshift,mass
    
    def TrackMainHaloProgenitorsOfFOF(self,TrackHaloID,SnapNum):
        index=np.where(self.SubhaloNr[np.where(self.SnapNum==SnapNum)[0]]==TrackHaloID)[0]
        print(index,TrackHaloID,self.SubhaloNr[index],self.GroupNr[index],self.SubhaloPos[index])
        # index=HaloID
        
        FOFindex=self.TreeFirstHaloInFOFgroup[index]   # This selects the main subhalo of the FOF group
        print(FOFindex)
        MainBranchIndex=FOFindex

        self.pos_main_branch=np.array([])
        self.vel_main_branch=np.array([])
        self.mass_main_branch=np.array([])        
        self.time_main_branch=np.array([])

        num_nodes=0
        while self.TreeMainProgenitor[MainBranchIndex]!=-1:
            self.pos_main_branch=np.append(self.pos_main_branch,self.SubhaloPos[MainBranchIndex]*self.Time[self.SnapNum[MainBranchIndex]])
            self.vel_main_branch=np.append(self.vel_main_branch,self.SubhaloVel[MainBranchIndex]*np.sqrt(self.Time[self.SnapNum[MainBranchIndex]]))
            self.mass_main_branch=np.append(self.mass_main_branch,self.GroupMass[MainBranchIndex])            
            self.time_main_branch=np.append(self.time_main_branch,self.Time[self.SnapNum[MainBranchIndex]])
            num_nodes+=1
            MainBranchIndex=self.TreeMainProgenitor[MainBranchIndex]
        print("Number of nodes of main branch: %d"%num_nodes)
        self.pos_main_branch=self.pos_main_branch.reshape(num_nodes,3)
        self.vel_main_branch=self.vel_main_branch.reshape(num_nodes,3)
        print("Earliest expansion factor: %f"%self.time_main_branch[num_nodes-1])

        # return time_main_branch,mass_main_branch,pos_main_branch,vel_main_branch

        num_in_FOF=0

        FOFindex=self.TreeNextHaloInFOFgroup[FOFindex] # This selects the first (i.e. not main) subhalo in the FOF group now
        
        subhalo_pos=[]
        subhalo_vel=[]
        subhalo_mass=[]
        subhalo_time=[]
        
        while self.TreeNextHaloInFOFgroup[FOFindex]!=-1:
            Treeindex=FOFindex

            dpos=np.array([])
            dvel=np.array([])
            mass=np.array([])
            time=np.array([])

            k=0
            
            while self.TreeMainProgenitor[Treeindex]!=-1:
                if np.any(self.Time[self.SnapNum[Treeindex]]==self.time_main_branch)==True:
                    j=np.where(self.time_main_branch==self.Time[self.SnapNum[Treeindex]])[0]
                    pos_main_subhalo=self.pos_main_branch[j]
                    pos_subhalo=self.SubhaloPos[Treeindex]*self.Time[self.SnapNum[Treeindex]]
                    d=pos_subhalo-pos_main_subhalo
                    d=np.where(d>0.5*self.BoxSize*self.Time[self.SnapNum[Treeindex]],d-self.BoxSize*self.Time[self.SnapNum[Treeindex]],d)
                    d=np.where(d<-0.5*self.BoxSize*self.Time[self.SnapNum[Treeindex]],d+self.BoxSize*self.Time[self.SnapNum[Treeindex]],d)
                    dpos=np.append(dpos,d)
                                                                    
                    vel_main_subhalo=self.vel_main_branch[j]
                    vel_subhalo=self.SubhaloVel[Treeindex]*np.sqrt(self.Time[self.SnapNum[Treeindex]])
                    dvel=np.append(dvel,vel_subhalo-vel_main_subhalo)
                    radius=np.sqrt(np.sum(dpos*dpos))
                    vrad=np.sum(dpos*dvel)/radius
                    mass=np.append(mass,self.SubhaloMass[Treeindex])
                    time=np.append(time,self.Time[self.SnapNum[Treeindex]])
                    k+=1
                Treeindex=self.TreeMainProgenitor[Treeindex]
            subhalo_time.append(time)
            subhalo_mass.append(mass)
            subhalo_pos.append(dpos.reshape(k,3))
            subhalo_vel.append(dvel.reshape(k,3))
 
            num_in_FOF+=1
            FOFindex=self.TreeNextHaloInFOFgroup[FOFindex]

        return num_in_FOF,subhalo_time,subhalo_mass,subhalo_pos,subhalo_vel
    
    def FindSubhaloIndexInMergerTree(self,ObjectID,ObjectType,SnapNum):
        if ObjectType=='Group':
            SubhaloID=self.GroupFirstSub[ObjectID]
        else:
            SubhaloID=ObjectID
        SubhaloIndex=np.where(self.SnapNum==SnapNum)[0][np.where(self.SubhaloNr[np.where(self.SnapNum==SnapNum)[0]]==SubhaloID)[0]][0]
        TreeID=self.TreeHalosID[SubhaloIndex]
        return SubhaloIndex,TreeID
                    
    def GetMergerTreeIDBounds(self,TreeID):
        TreeIDStart=tree.TreeOffset[TreeID]
        TreeIDEnd=TreeIDStart+tree.TreeLength[TreeID]
        return TreeIDStart,TreeIDEnd        

    """ Trees for the subhalos within a FOF group are grouped together in the same tree structure.
    """
    def GetMergerTree(self,TreeID):
        TreeIDStart,TreeIDFinish=self.GetMergerTreeIDBounds(TreeID)
        tree_index=self.TreeHalosIndex[TreeIDStart:TreeIDFinish]
        tree_prog=self.TreeProgenitor[TreeIDStart:TreeIDFinish]
        tree_des=self.TreeDescendant[TreeIDStart:TreeIDFinish]
        tree_nextdes=self.TreeNextDescendant[TreeIDStart:TreeIDFinish]
        tree_nextprog=self.TreeNextProgenitor[TreeIDStart:TreeIDFinish]
        tree_grpnum=self.GrpNr[TreeIDStart:TreeIDFinish]
        tree_subnum=self.SubhaloNr[TreeIDStart:TreeIDFinish]
        tree_snapnum=self.SnapNum[TreeIDStart:TreeIDFinish]
        
        edge_list=[]
        snap_list=[]
        coords_list={}
        
        index=0
        i=0
        coords_list[index]=[i,tree_snapnum[index]]

        while tree_prog[index]!=-1:
#            coords_list[index]=[i,tree_snapnum[index]]
            edge_list+=[[index,tree_prog[index]]]
            snap_list+=[[index,tree_prog[index],tree_snapnum[index],tree_snapnum[tree_prog[index]]]]
            num_progenitors=0
            prindex=tree_prog[index]
            j=i+0.1
            while tree_nextprog[prindex]!=-1:
                edge_list+=[[index,tree_nextprog[prindex]]]        
                snap_list+=[[index,tree_nextprog[prindex],tree_snapnum[index],tree_snapnum[tree_nextprog[prindex]]]]                
                num_progenitors+=1
                prindex=tree_nextprog[prindex]         
                if prindex!=-1:
                    coords_list[prindex]=[j,tree_snapnum[prindex]]
                j+=0.1
                    
#            print(tree.Time[tree_snapnum[index]],index,tree_subnum[index],tree_subnum[tree_des[index]],tree_nextprog[index],tree_subnum[tree_des[tree_nextprog[index]]],num_progenitors)
            index=tree_prog[index]
#            i+=1
            
            if index!=-1:
                coords_list[index]=[i,tree_snapnum[index]]
            
            if tree_snapnum[index]<55:
                break

        
        return edge_list,coords_list,snap_list
    
    def plot_orbits(self,HaloID,f_size=1,amin=0.2,amax=1.0,rmin=0.0,rmax=2.0,plot_redshift=False,**kwargs):
        self.ReadMergerTreeFile()
        n,a,m,pos,vel=self.TrackMainHaloProgenitorsOfFOF(HaloID,SnapNum)
        rhocrit_0=27.755
        delta_vir=200
        rhocrit_z=rhocrit_0*np.sqrt(self.Omega0*self.time_main_branch**-3.+self.OmegaLambda)
        r200_z=f_size*(3*self.mass_main_branch/4./np.pi/delta_vir/rhocrit_z)**(1./3.)
        r200_z*=self.time_main_branch

        fig=plt.figure(figsize=(8,8))
        ax=fig.add_subplot(111)

        if plot_redshift==True:                
            ax.set_xlim(1./amin-1,1./amax-1)
        else:
            ax.set_xlim(amin,amax)
        ax.set_ylim(rmin,rmax)

        for i in range(n):
            r=np.array([])
            for j in range(len(a[i])):
                r=np.append(r,np.sqrt(np.sum(pos[i][j]**2)))
            if plot_redshift==True:
                ax.plot(1./a[i]-1,r)
            else:
                ax.plot(a[i],r)
        if plot_redshift==True:                
            ax.plot(1./self.time_main_branch-1,r200_z,lw=3,ls=":",color="black")
            ax.set_xlabel("Redshift")            
        else:
            ax.plot(self.time_main_branch,r200_z,lw=3,ls=":",color="black")
            ax.set_xlabel("Expansion Factor")
            
        ax.set_ylabel("Subhalo Radius [pMpc/h]")
        return fig
        
    def plot_halo(self,SubhaloID,SnapNum,f_size=1,**kwargs):
        self.ReadMergerTreeFile()
        self.ReadHaloCatalogue()
        self.ReadSnapshot()
        n,a,m,pos,vel=self.TrackMainHaloProgenitorsOfFOF(SubhaloID,SnapNum)

        GroupID=self.SubhaloGroupNr[SubhaloID]
        print(n,self.GroupNsubs[GroupID])        
        pos_subhalos=np.array([])
        num_subs=0
        for i in range(len(pos)):
            if len(pos[i])>0:
                num_subs+=1
                j=np.where(np.abs(a[i]-self.GroupAscale[GroupID])<1.e-4)[0]
                pos_subhalos=np.append(pos_subhalos,pos[i][j])
        pos_subhalos=pos_subhalos.reshape(num_subs,3)
        
        r200_z=self.GroupR200[GroupID]*self.GroupAscale[GroupID]

        virial_circle=plt.Circle((0,0),self.GroupR200[GroupID],color="black",fill=False,lw=3,ls=":")

        fig=plt.figure(figsize=(8,8))        
        ax=fig.add_subplot(111)
        ax.set_xlim(f_size*r200_z,-f_size*r200_z)
        ax.set_ylim(f_size*r200_z,-f_size*r200_z)
        
        ax.scatter((self.pos-self.GroupPos[GroupID])[:,0],(self.pos-self.GroupPos[GroupID])[:,1],s=0.01,color="grey",alpha=0.3)        

        # istart=self.SubhaloOffsetType[SubaloID][1]
        # ifinish=istart+np.sum(tree.SubhaloLen[np.where(tree.SubhaloGroupNr==tree.SubhaloGroupNr[HaloID])[0]])
        
        # ax.scatter(self.pos[istart:ifinish,1]-self.SubPos[HaloID][1],self.pos[istart:ifinish,2]-self.SubPos[HaloID][2],s=0.01,color="blue")        
        ax.scatter(pos_subhalos[:,1],pos_subhalos[:,2],marker="*",color="red")
        # ax.scatter(self.SubPos[0:self.GroupNsubs[0]-1,1]-self.GroupPos[HaloID][1],self.SubPos[0:self.GroupNsubs[0]-1,2]-self.GroupPos[HaloID][2],s=1,color="green")

        plt.gca().add_patch(virial_circle)

        ax.set_xlabel("X [Mpc/h]")
        ax.set_ylabel("Y [Mpc/h]")
        return fig
        
# def GetMainBranchLength(SnapShotNr):
#     ids=np.where(tree.SnapNum==SnapShotNr)[0]
#     ids=np.where(tree.TreeHalosIndex[ids]==0)[0]
#     num_halos=len(ids)
#     length_branch=np.array([],dtype=np.int32)
#     for i in range(num_halos-1):
#         index=ids[i]
#         j=0
#         while tree.TreeMainProgenitor[index] != -1:
#             j=j+1 
#             index=tree.TreeMainProgenitor[index]
#         length_branch=np.append(length_branch,j)
#     return length_branch
    
