import h5py
import numpy as np
import os
from typing import Union, List, Optional

class read_hdf5:
    def __init__(self):
        pass
    def read_hdf5_snapshot(self,snapfilename: Optional[str] = None,convention: Optional[str] = None):
        """Main function to read HDF5 snapshot data."""
        # Determine filename
        filename = self._determine_filename(snapfilename)
        print('Reading data from %s' % filename)
        
        if convention is not None:
            self.convention = convention
        
        # Read header and setup
        with h5py.File(filename, 'r') as f:
            NameOfMassBlock = self._read_header_info(f)
            self._read_cosmological_params(f)
            self._check_potential_availability(f)
        
        # Calculate total particles and print summary
        NumPart = np.sum(self.NumPart_Total)
        self._print_summary(NumPart)
        
        # Allocate memory
        self._allocate_memory(NumPart)
        extra_flags = self._allocate_extra_block_memory()
        
        # Calculate offsets
        istart, ifinish = self._calculate_offsets()
        
        # Read particle data
        if self.NumFiles > 1:
            self._read_particle_data_multiple_files(self.snapfilename, NameOfMassBlock,
                                                   istart, ifinish, extra_flags)
        else:
            self._read_particle_data_single_file(filename, NameOfMassBlock,
                                                istart, ifinish, extra_flags)


    def _determine_filename(self,snapfilename):
        """Determine the correct filename for the snapshot."""
        filename = snapfilename + '.hdf5'
        if not os.path.exists(filename):
            filename = snapfilename + '.0.hdf5'
        return filename

    def _read_header_info(self, f):
        """Extract header information from the HDF5 file."""
        self.NumFiles = f['Header'].attrs['NumFilesPerSnapshot']
        self.NumPart_Total = f['Header'].attrs['NumPart_Total'][()]
        self.NumPartType = np.shape(self.NumPart_Total)[0]
        self.MassTable = f['Header'].attrs['MassTable'][()]
        
        # Determine mass block name
        NameOfMassBlock = 'Mass'
        ii = np.where(f['Header'].attrs['MassTable'][()] == 0)[0]
        for indx in ii:
            if f['Header'].attrs['NumPart_Total'][indx] > 0:
                if 'Masses' in list(f['PartType%d' % indx].keys()):
                    NameOfMassBlock = 'Masses'
        
        return NameOfMassBlock

    
    def _read_cosmological_params(self, f):
        """Read cosmological parameters based on simulation convention."""
        if self.convention == 'SWIFT':
            self.ScaleFactor = f['Header'].attrs['Scale-factor']
            self.BoxSize = f['Header'].attrs['BoxSize'][()][0]
            self.OmegaDM = f['Cosmology'].attrs['Omega_cdm']
            self.OmegaBar = f['Cosmology'].attrs['Omega_b']
            self.Omega0 = self.OmegaDM + self.OmegaBar
            self.OmegaLambda = f['Cosmology'].attrs['Omega_lambda']
            self.HubbleParam = f['Cosmology'].attrs['h']
        elif self.convention in ['GADGET4', 'AREPO']:
            self.ScaleFactor = f['Header'].attrs['Time']
            self.BoxSize = f['Header'].attrs['BoxSize']
            self.Omega0 = f['Parameters'].attrs['Omega0']
            self.OmegaLambda = f['Parameters'].attrs['OmegaLambda']
            self.HubbleParam = f['Parameters'].attrs['HubbleParam']
        else:
            self.ScaleFactor = f['Header'].attrs['Time']
            self.BoxSize = f['Header'].attrs['BoxSize']
            self.Omega0 = f['Header'].attrs['Omega0']
            self.OmegaLambda = f['Header'].attrs['OmegaLambda']
            self.HubbleParam = f['Header'].attrs['HubbleParam']

    def _write_cosmological_params(self, f):
        """Read cosmological parameters based on simulation convention."""
        if self.convention == 'SWIFT':
            f['Header'].attrs['Scale-factor'] = self.ScaleFactor
            f['Header'].attrs['BoxSize'][()][0] = self.BoxSize
            f['Cosmology'].attrs['Omega_cdm'] = self.Omega0 - self.OmegaBar
            f['Cosmology'].attrs['Omega_b'] = self.OmegaBar
            f['Cosmology'].attrs['Omega_lambda'] = self.OmegaLambda
            f['Cosmology'].attrs['h'] = self.HubbleParam
        elif self.convention in ['GADGET4', 'AREPO']:
            self.ScaleFactor = f['Header'].attrs['Time'] = self.ScaleFactor
            self.BoxSize = f['Header'].attrs['BoxSize'] = self.BoxSize
            self.Omega0 = f['Parameters'].attrs['Omega0'] = self.Omega0
            self.OmegaLambda = f['Parameters'].attrs['OmegaLambda'] = self.OmegaLambda
            self.HubbleParam = f['Parameters'].attrs['HubbleParam'] = self.HubbleParam
        else:
            f['Header'].attrs['Time'] = self.ScaleFactor
            f['Header'].attrs['BoxSize'] = self.BoxSize
            f['Header'].attrs['Omega0'] = self.Omega0
            f['Header'].attrs['OmegaLambda'] = self.OmegaLambda
            f['Header'].attrs['HubbleParam'] = self.HubbleParam

    def _check_potential_availability(self, f):
        """Check if potential data is available."""
        self.ispotential = False
        if 'PartType1' in list(f.keys()):
            if 'Potential' in f['PartType1'].keys():
                self.ispotential = True

    def _print_summary(self, NumPart):
        """Print summary information about the snapshot."""
        print('Simulation scale factor: %lf' % self.ScaleFactor)
        if self.NumFiles > 1:
            print('Data is split across %d files' % self.NumFiles)
        
        print('Number of particles: %010d' % NumPart)
        print('Number of particle types: %d' % self.NumPartType)
        
        idx_with_mass = np.where(self.MassTable > 0)[0]
        NumPart_InMassBlock = np.sum(self.NumPart_Total[idx_with_mass])
        if NumPart_InMassBlock > 0:
            print('Number of particles in mass block: %010d' % NumPart_InMassBlock)
        
        if getattr(self,'hires_only',False):
            hires_particles = (self.NumPart_Total[self.gas_type] +
                             self.NumPart_Total[self.dm_type] +
                             self.NumPart_Total[self.star_type] +
                             self.NumPart_Total[self.bh_type])
            print('Number of HIRES particles: %010d' % hires_particles)

    def _allocate_memory(self, NumPart):
        """Allocate memory for particle data arrays."""
        # Position arrays
        if self.positions_type == 'float64':
            self.pos = np.ndarray(shape=(NumPart, 3), dtype=np.float64)
        else:
            self.pos = np.ndarray(shape=(NumPart, 3), dtype=np.float32)
        
        # Basic particle data
        if not self.positions_only:
            self.vel = np.ndarray(shape=(NumPart, 3))
            if self.pids_type == 32:
                self.pids = np.ndarray(shape=(NumPart), dtype=np.uint32)
            else:
                self.pids = np.ndarray(shape=(NumPart), dtype=np.uint64)
            self.mass = np.ndarray(shape=(NumPart))
        
        # Gas particle specific arrays
        if self.NumPart_Total[0] > 0 and not self.positions_only:
            self.u = np.ndarray(shape=(self.NumPart_Total[0]))
            self.rho = np.ndarray(shape=(self.NumPart_Total[0]))
            if self.convention != 'Arepo':
                self.smoothinglength = np.ndarray(shape=(self.NumPart_Total[0]))
        
        # Particle type array
        if self.get_ptypes:
            self.ptype = np.ones(shape=(NumPart), dtype=np.int32)
            ioffset = np.zeros(shape=(self.NumPartType + 1), dtype=np.uint64)
            ioffset[1:] = np.cumsum(self.NumPart_Total)
            for i in range(self.NumPartType):
                self.ptype[ioffset[i]:ioffset[i+1]] = self.ptype[ioffset[i]:ioffset[i+1]] * i
        
        # Potential array
        if self.ispotential:
            self.potential = np.ndarray(shape=(NumPart))

    def _allocate_extra_block_memory(self):
        """Allocate memory for extra data blocks."""
        extra_flags = {}
        
        if len(self.extra_blocks) > 0:
            print('Loading extra blocks: %s' % self.extra_blocks)
            
            if 'AGE' in self.extra_blocks:
                self.stellarage = np.ndarray(shape=(self.NumPart_Total[self.star_type]))
                extra_flags['isstellarage'] = True
            
            if 'Z' in self.extra_blocks:
                self.gas_metallicity = np.ndarray(shape=(self.NumPart_Total[self.gas_type]))
                self.stellar_metallicity = np.ndarray(shape=(self.NumPart_Total[self.star_type]))
                extra_flags['ismetallicity'] = True
            
            if 'SFR' in self.extra_blocks:
                self.gas_sfr = np.ndarray(shape=(self.NumPart_Total[self.gas_type]))
                extra_flags['issfr'] = True
            
            if 'STELLARGENS' in self.extra_blocks:
                self.stellargen = np.ndarray(shape=(self.NumPart_Total[self.star_type]), dtype=np.int32)
                extra_flags['isstellargens'] = True
            
            if 'INIT_MASS' in self.extra_blocks:
                self.stellarinitmass = np.ndarray(shape=(self.NumPart_Total[self.star_type]), dtype=np.float32)
                extra_flags['isstellarinitmass'] = True
            
            if 'POT' in self.extra_blocks:
                self.potential = np.ndarray(shape=(np.sum(self.NumPart_Total)))
                extra_flags['ispotential'] = True
        
        return extra_flags

    def _calculate_offsets(self):
        """Calculate starting offsets for each particle type."""
        istart = np.zeros(self.NumPartType, dtype=np.uint64)
        offset = 0
        
        for i in range(self.NumPartType):
            if self.NumPart_Total[i] > 0:
                istart[i] = offset
            offset += self.NumPart_Total[i]
        
        ifinish = np.copy(istart)
        print('istart:', istart)
        
        return istart, ifinish

    def _read_particle_data_single_file(self, filename, NameOfMassBlock, istart, ifinish, extra_flags):
        """Read particle data from a single HDF5 file."""
        with h5py.File(filename, 'r') as f:
            jstart = 0
            
            for itype in range(self.NumPartType):
                if self.hires_only and itype in self.not_hires_ptypes:
                    print('Skipping Particle Type %d' % itype)
                    continue
                
                if self.NumPart_Total[itype] > 0:
                    ifinish[itype] = istart[itype] + self.NumPart_Total[itype]
                    
                    # Read basic particle data
                    self._read_basic_particle_data(f, itype, istart[itype], ifinish[itype], NameOfMassBlock)
                    
                    # Read type-specific data
                    if itype == self.gas_type:
                        jstart = self._read_gas_data(f, itype, istart[itype], ifinish[itype],
                                                   jstart, extra_flags)
                    elif itype == self.star_type:
                        jstart = self._read_star_data(f, itype, istart[itype], ifinish[itype],
                                                    jstart, extra_flags)
                    
                    istart[itype] = ifinish[itype]

    def _read_particle_data_multiple_files(self, base_filename, NameOfMassBlock, istart, ifinish, extra_flags):
        """Read particle data from multiple HDF5 files."""
        for i in range(self.NumFiles):
            filename = base_filename + '.%d.hdf5' % i
            print('Reading in file %s...' % filename)
            
            with h5py.File(filename, 'r') as f:
                NumPart_ThisFile = f['Header'].attrs['NumPart_ThisFile'][()]
                jstart = 0
                
                for itype in range(self.NumPartType):
                    if self.hires_only and itype in self.not_hires_ptypes:
                        print('Skipping Particle Type %d' % itype)
                        continue
                    
                    if NumPart_ThisFile[itype] > 0:
                        ifinish[itype] = istart[itype] + NumPart_ThisFile[itype]
                        
                        # Read basic particle data
                        self._read_basic_particle_data(f, itype, istart[itype], ifinish[itype], NameOfMassBlock)
                        
                        # Read type-specific data
                        if itype == self.gas_type:
                            jstart = self._read_gas_data(f, itype, istart[itype], ifinish[itype],
                                                       jstart, extra_flags)
                        elif itype == self.star_type:
                            jstart = self._read_star_data(f, itype, istart[itype], ifinish[itype],
                                                        jstart, extra_flags, NumPart_ThisFile[itype])
                        
                        istart[itype] = ifinish[itype]

    def _read_basic_particle_data(self, f, itype, istart, ifinish, NameOfMassBlock):
        """Read basic particle data (positions, velocities, IDs, masses, potential)."""
        # Positions
        self.pos[istart:ifinish] = f['PartType%d/Coordinates' % itype][()]
        
        if self.positions_only:
            return
        
        # Velocities and IDs
        self.vel[istart:ifinish] = f['PartType%d/Velocities' % itype][()]
        self.pids[istart:ifinish] = f['PartType%d/ParticleIDs' % itype][()]
        
        # Masses
        if self.MassTable[itype] == 0:
            if NameOfMassBlock == 'Mass':
                self.mass[istart:ifinish] = f['PartType%d/Mass' % itype][()]
            else:
                if itype == 5 and self.convention == 'SWIFT':
                    self.mass[istart:ifinish] = f['PartType%d/DynamicalMasses' % itype][()]
                else:
                    self.mass[istart:ifinish] = f['PartType%d/Masses' % itype][()]
        else:
            self.mass[istart:ifinish] = self.MassTable[itype]
        
        # Potential
        if self.ispotential:
            self.potential[istart:ifinish] = f['PartType%d/Potential' % itype][()]

    def _read_gas_data(self, f, itype, istart, ifinish, jstart, extra_flags):
        """Read gas particle specific data."""
        if self.convention == 'SWIFT':
            self.u[istart:ifinish] = f['PartType%d/InternalEnergies' % itype][()]
            self.rho[istart:ifinish] = f['PartType%d/Densities' % itype][()]
        else:
            self.u[istart:ifinish] = f['PartType%d/InternalEnergy' % itype][()]
            self.rho[istart:ifinish] = f['PartType%d/Density' % itype][()]
            if self.convention != 'Arepo':
                self.smoothinglength[istart:ifinish] = f['PartType%d/SmoothingLength' % itype][()]
        
        # Extra gas data
        if extra_flags.get('ismetallicity', False):
            if f['PartType%d/Metallicity' % itype].ndim > 1:
                self.gas_metallicity[istart:ifinish] = f['PartType%d/Metallicity' % itype][:, 0]
            else:
                self.gas_metallicity[istart:ifinish] = f['PartType%d/Metallicity' % itype][()]
        
        if extra_flags.get('issfr', False):
            self.gas_sfr[istart:ifinish] = f['PartType%d/StarFormationRate' % itype][()]
        
        return jstart

    def _read_star_data(self, f, itype, istart, ifinish, jstart, extra_flags, num_particles=None):
        """Read star particle specific data."""
        if num_particles is None:
            jfinish = jstart + self.NumPart_Total[itype]
        else:
            jfinish = jstart + num_particles
        
        if extra_flags.get('ismetallicity', False):
            if f['PartType%d/Metallicity' % itype].ndim > 1:
                self.stellar_metallicity[jstart:jfinish] = f['PartType%d/Metallicity' % itype][:, 0]
            else:
                self.stellar_metallicity[jstart:jfinish] = f['PartType%d/Metallicity' % itype][()]
        
        if extra_flags.get('isstellarage', False):
            self.stellarage[jstart:jfinish] = f['PartType%d/StellarFormationTime' % itype][()]
        
        if extra_flags.get('isstellargens', False):
            self.stellargen[jstart:jfinish] = f['PartType%d/ID_Generations' % itype][()]
        
        if extra_flags.get('isstellarinitmass', False):
            self.stellarinitmass[jstart:jfinish] = f['PartType%d/StellarInitMass' % itype][()]
        
        return jfinish

