"""
Module: snapio_hdf5 contains ReadHDF5 and WriteHDF5, which provide tools to read and write HDF5 files given
a particular convention.
"""
import h5py
import numpy as np
import os
from typing import Union, List, Optional

class read_hdf5:
    """
    ReadHDF5 - a class to read HDF5 files in GADGET, Arepo, and SWIFT format.
    Usage:
        from snapio_hdf5 import read_hdf5
        reader = ReadHDF5()
        self._transfer_attributes_to_reader(reader)
        reader.read_hdf5_snapshot(snapfilename=filename,convention=convention)
        self._transfer_attributes_from_reader(reader)
    """

    def __init__(self):
        pass
        
    def read_hdf5_snapshot(self,snapfilename: Optional[str] = None,convention: Optional[str] = None):
        """Main function to read HDF5 snapshot data."""
        # Determine filename
        filename = self._determine_filename(snapfilename)
        print(f"Reading data from {filename}")
        
        if convention is not None:
            self.convention = convention
        # Read header and setup
        with h5py.File(filename, 'r') as f:
            name_of_mass_block = self._read_header_info(f)
            self._read_cosmological_params(f)
            self._check_potential_availability(f)
        
        # Calculate total particles and print summary
        num_part = np.sum(self.num_part_total)
        self._print_summary(num_part)
        
        # Allocate memory
        self._allocate_memory(num_part)
        extra_flags = self._allocate_extra_block_memory()
        
        # Calculate offsets
        istart, ifinish = self._calculate_offsets()
        
        # Read particle data
        if self.num_files > 1:
            self._read_particle_data_multiple_files(snapfilename, name_of_mass_block,
                                                   istart, ifinish, extra_flags)
        else:
            self._read_particle_data_single_file(filename, name_of_mass_block,
                                                istart, ifinish, extra_flags)


    def _determine_filename(self,snapfilename):
        """Determine the correct filename for the snapshot."""
        if os.path.exists(snapfilename):
            return snapfilename
        filename = snapfilename + '.hdf5'
        if not os.path.exists(filename):
            filename = snapfilename + '.0.hdf5'
        return filename

    def _read_header_info(self, f):
        """Extract header information from the HDF5 file."""
        self.num_files = f['Header'].attrs['NumFilesPerSnapshot']
        self.num_part_total = f['Header'].attrs['NumPart_Total'][()]
        self.num_type = self.num_part_total.size
        self.mass_table = f['Header'].attrs['MassTable'][()]
        
        # Determine mass block name
        name_of_mass_block = 'Mass'
        ii = np.where(self.mass_table == 0)[0]
        for indx in ii:
            if self.num_part_total[indx] > 0:
                if 'Masses' in list(f['PartType%d' % indx].keys()):
                    name_of_mass_block = 'Masses'
        
        return name_of_mass_block

    def _read_cosmological_params(self, f):
        """Read cosmological parameters based on simulation convention."""
        if self.convention == 'SWIFT':
            self.scale_factor = f['Header'].attrs['Scale-factor']
            self.dimension = f['Header'].attrs['Dimension']
            self.box_size = f['Header'].attrs['BoxSize'][()][0]
            self.omega_dm = f['Cosmology'].attrs['Omega_cdm']
            self.omega_bar = f['Cosmology'].attrs['Omega_b']
            self.omega_0 = self.omega_dm + self.omega_bar
            self.omega_lambda = f['Cosmology'].attrs['Omega_lambda']
            self.hubble_param = f['Cosmology'].attrs['h']
        elif self.convention in ['GADGET4', 'AREPO']:
            self.scale_factor = f['Header'].attrs['Time']
            self.box_size = f['Header'].attrs['BoxSize']
            self.omega_0 = f['Parameters'].attrs['Omega0']
            self.omega_lambda = f['Parameters'].attrs['OmegaLambda']
            self.hubble_param = f['Parameters'].attrs['HubbleParam']
        else:
            self.scale_factor = f['Header'].attrs['Time']
            self.box_size = f['Header'].attrs['BoxSize']
            self.omega_0 = f['Header'].attrs['Omega0']
            self.omega_lambda = f['Header'].attrs['OmegaLambda']
            self.hubble_param = f['Header'].attrs['HubbleParam']

    def _check_potential_availability(self, f):
        """Check if potential data is available."""
        self.ispotential = False
        if 'PartType1' in list(f.keys()):
            if 'Potential' in f['PartType1'].keys():
                self.ispotential = True

    def _print_summary(self, num_part):
        """Print summary information about the snapshot."""
        print(f"Simulation scale factor: {self.scale_factor}")
        if self.num_files > 1:
            print(f"Data is split across {self.num_files} files")
        
        print(f"Number of particles: {num_part}")
        print(f"Number of particle types: {self.num_part_total.size}")
        
        idx_with_mass = np.where(self.mass_table > 0)[0]
        num_part_in_mass_block = np.sum(self.num_part_total[idx_with_mass])
        if num_part_in_mass_block > 0:
            print(f"Number of particles in mass block: {num_part_in_mass_block}")
        
        if getattr(self,'hires_only',False):
            hires_particles = (self.num_part_total[self.gas_type] +
                             self.num_part_total[self.dm_type] +
                             self.num_part_total[self.star_type] +
                             self.num_part_total[self.bh_type])
            print(f"Number of HIRES particles: {hires_particles}")

    def _allocate_memory(self, num_part):
        """Allocate memory for particle data arrays."""
        # Position arrays
        if self.positions_type == 'float64':
            self.pos = np.ndarray(shape=(num_part, 3), dtype=np.float64)
        else:
            self.pos = np.ndarray(shape=(num_part, 3), dtype=np.float32)
        
        # Basic particle data
        if not self.positions_only:
            self.vel = np.ndarray(shape=(num_part, 3))
            if self.pids_type == 32:
                self.pids = np.ndarray(shape=(num_part), dtype=np.uint32)
            else:
                self.pids = np.ndarray(shape=(num_part), dtype=np.uint64)
            self.mass = np.ndarray(shape=(num_part))
        
        # Gas particle specific arrays
        if self.num_part_total[0] > 0 and not self.positions_only:
            self.u = np.ndarray(shape=(self.num_part_total[0]))
            self.rho = np.ndarray(shape=(self.num_part_total[0]))
            self.smoothinglength = np.ndarray(shape=(self.num_part_total[0]))
        
        # Particle type array
        if self.get_ptypes:
            self.ptype = np.ones(shape=(num_part), dtype=np.int32)
            ioffset = np.zeros(shape=(self.num_type + 1), dtype=np.uint64)
            ioffset[1:] = np.cumsum(self.num_part_total)
            for i in range(self.num_type):
                self.ptype[ioffset[i]:ioffset[i+1]] = self.ptype[ioffset[i]:ioffset[i+1]] * i
        
        # Potential array
        if self.ispotential:
            self.potential = np.ndarray(shape=(num_part))

    def _allocate_extra_block_memory(self):
        """Allocate memory for extra data blocks."""
        extra_flags = {}
        
        if len(self.extra_blocks) > 0:
            print(f"Loading extra blocks: {self.extra_blocks}")
            if 'AGE' in self.extra_blocks:
                self.stellarage = np.ndarray(shape=(self.num_part_total[self.star_type]))
                extra_flags['isstellarage'] = True
            if 'Z' in self.extra_blocks:
                self.gas_metallicity = np.ndarray(shape=(self.num_part_total[self.gas_type]))
                self.stellar_metallicity = np.ndarray(shape=(self.num_part_total[self.star_type]))
                extra_flags['ismetallicity'] = True
            if 'SFR' in self.extra_blocks:
                self.gas_sfr = np.ndarray(shape=(self.num_part_total[self.gas_type]))
                extra_flags['issfr'] = True
            if 'STELLARGENS' in self.extra_blocks:
                self.stellargen = np.ndarray(shape=(self.num_part_total[self.star_type]), dtype=np.int32)
                extra_flags['isstellargens'] = True
            if 'INIT_MASS' in self.extra_blocks:
                self.stellarinitmass = np.ndarray(shape=(self.num_part_total[self.star_type]), dtype=np.float32)
                extra_flags['isstellarinitmass'] = True
            if 'POT' in self.extra_blocks:
                self.potential = np.ndarray(shape=(np.sum(self.num_part_total)))
                extra_flags['ispotential'] = True
                self.ispotential = True
        
        return extra_flags

    def _calculate_offsets(self):
        """Calculate starting offsets for each particle type."""
        istart = np.zeros(self.num_type, dtype=np.uint64)
        offset = 0
        for i in range(self.num_type):
            if self.num_part_total[i] > 0:
                istart[i] = offset
            offset += self.num_part_total[i]
        ifinish = np.copy(istart)
        return istart, ifinish

    def _read_particle_data_single_file(self, filename, name_of_mass_block, istart, ifinish, extra_flags):
        """Read particle data from a single HDF5 file."""
        with h5py.File(filename, 'r') as f:
            jstart = 0
            for itype in range(self.num_type):
                if self.hires_only and itype in self.not_hires_ptypes:
                    print(f"Skipping Particle Type {itype}")
                    continue
                if self.num_part_total[itype] > 0:
                    ifinish[itype] = istart[itype] + self.num_part_total[itype]
                    # Read basic particle data
                    self._read_basic_particle_data(f, itype, istart[itype], ifinish[itype], name_of_mass_block)
                    # Read type-specific data
                    if itype == self.gas_type:
                        jstart = self._read_gas_data(f, itype, istart[itype], ifinish[itype],
                                                   jstart, extra_flags)
                    elif itype == self.star_type:
                        jstart = self._read_star_data(f, itype, istart[itype], ifinish[itype],
                                                    jstart, extra_flags)
                    istart[itype] = ifinish[itype]

    def _read_particle_data_multiple_files(self, base_filename, name_of_mass_block, istart, ifinish, extra_flags):
        """Read particle data from multiple HDF5 files."""
        for i in range(self.num_files):
            filename = base_filename + '.%d.hdf5' % i
            print(f"Reading in file {filename}..")
            with h5py.File(filename, 'r') as f:
                num_part_this_file = f['Header'].attrs['num_part_this_file'][()]
                jstart = 0
                for itype in range(self.num_type):
                    if self.hires_only and itype in self.not_hires_ptypes:
                        print('Skipping Particle Type %d' % itype)
                        continue
                    if num_part_this_file[itype] > 0:
                        ifinish[itype] = istart[itype] + num_part_this_file[itype]
                        # Read basic particle data
                        self._read_basic_particle_data(f, itype, istart[itype], ifinish[itype], name_of_mass_block)
                        # Read type-specific data
                        if itype == self.gas_type:
                            jstart = self._read_gas_data(f, itype, istart[itype], ifinish[itype],
                                                    jstart, extra_flags)
                        elif itype == self.star_type:
                            jstart = self._read_star_data(f, itype, istart[itype], ifinish[itype],
                                                    jstart, extra_flags, num_part_this_file[itype])
                        istart[itype] = ifinish[itype]

    def _read_basic_particle_data(self, f, itype, istart, ifinish, name_of_mass_block):
        """Read basic particle data (positions, velocities, IDs, masses, potential)."""
        # Positions
        self.pos[istart:ifinish] = f['PartType%d/Coordinates' % itype][()]
        if self.positions_only:
            return
        # Velocities and IDs
        self.vel[istart:ifinish] = f['PartType%d/Velocities' % itype][()]
        self.pids[istart:ifinish] = f['PartType%d/ParticleIDs' % itype][()]
        # Masses
        if self.mass_table[itype] == 0:
            if name_of_mass_block == 'Mass':
                self.mass[istart:ifinish] = f['PartType%d/Mass' % itype][()]
            else:
                if itype == 5 and self.convention == 'SWIFT':
                    self.mass[istart:ifinish] = f['PartType%d/DynamicalMasses' % itype][()]
                else:
                    self.mass[istart:ifinish] = f['PartType%d/Masses' % itype][()]
        else:
            self.mass[istart:ifinish] = self.mass_table[itype]
        
        # Potential
        if self.ispotential:
            g = f[f"PartType{itype}"]
            for key in ("Potential", "Potentials"):
                if key in g:
                    print(f"Reading Potential for Particle Type {itype} as {key}")
                    self.potential[istart:ifinish] = g[key][()]
                    break
        else:
            raise KeyError("No Potential field found")

    def _read_gas_data(self, f, itype, istart, ifinish, jstart, extra_flags):
        """Read gas particle specific data."""
        ptype_group = f[f'PartType{itype}']
        field_name = 'InternalEnergy' if 'InternalEnergy' in ptype_group else 'InternalEnergies'
        self.u[istart:ifinish] = ptype_group[field_name][()]

        if not self.isics:
            field_name = 'Density' if 'Density' in ptype_group else 'Densities'
            self.rho[istart:ifinish] = ptype_group[field_name][()]

        if self.convention.lower() != 'arepo':
            self.smoothinglength[istart:ifinish] = f['PartType%d/SmoothingLengths' % itype][()]
        
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
            jfinish = jstart + self.num_part_total[itype]
        else:
            jfinish = jstart + num_particles
        
        if extra_flags.get('ismetallicity', False):
            if f['part_type%d/Metallicity' % itype].ndim > 1:
                self.stellar_metallicity[jstart:jfinish] = f['part_type%d/Metallicity' % itype][:, 0]
            else:
                self.stellar_metallicity[jstart:jfinish] = f['part_type%d/Metallicity' % itype][()]
        
        if extra_flags.get('isstellarage', False):
            self.stellarage[jstart:jfinish] = f['part_type%d/StellarFormationTime' % itype][()]
        
        if extra_flags.get('isstellargens', False):
            self.stellargen[jstart:jfinish] = f['part_type%d/ID_Generations' % itype][()]
        
        if extra_flags.get('isstellarinitmass', False):
            self.stellarinitmass[jstart:jfinish] = f['part_type%d/StellarInitMass' % itype][()]
        
        return jfinish

class write_hdf5:
    """
    write_hdf5 - a class to read HDF5 files in GADGET, Arepo, and SWIFT format.
    Usage:
        from snapio_hdf5 import write_hdf5
        writer = write_hdf5()
        self._transfer_attributes_to_writer(writer)
        writer.write_hdf5_snapshot(snapfilename=filename,convention=convention)
    """
    def __init__(self, output_convention, idx, idx_type):
        self.output_convention=output_convention
        self.idx = idx
        self.idx_type = idx_type
        if not hasattr(self, "name_of_u_block"):
            self.name_of_u_block = "InternalEnergy"
        if not hasattr(self, "name_of_mass_block"):
            self.name_of_mass_block = "Masses"
                
    def write_hdf5_snapshot(self, output_file):
        filename = output_file + ".hdf5"
        print(f"Writing data to {filename} in {self.output_convention.upper()} format.")
    
        with h5py.File(filename, "w") as f:
            self._write_header(f)
            self._write_particles(f,name_of_u_block=self.name_of_u_block,name_of_mass_block=self.name_of_mass_block)

    def _write_header(self, f) -> None:
        header = f.create_group("Header")
        base_attrs = {
            "NumFilesPerSnapshot": getattr(self,'num_files',1),
            "NumPart_Total": self.num_part_total,
            "NumPart_Total_HighWord": np.zeros(6, dtype=np.int32),
            "Flag_Entropy_ICs": 1,
            "MassTable": self.mass_table,
        }
        header.attrs.update(base_attrs)
        
        if self.output_convention.upper() == 'SWIFT':
            header.attrs.update({
                "Scale-factor": getattr(self,'scale_factor',1.),
                "Time": getattr(self,'time',0.),
                "BoxSize": np.full(3, getattr(self,'box_size',0.)),
                "Dimension": getattr(self,'dimension',3),
            })
        elif self.output_convention.upper() in ['GADGET4', 'AREPO']:
            redshift = 0.0 if self.scale_factor <= 0 else 1. / self.scale_factor - 1.
            header.attrs.update({
                "Time": getattr(self,'scale_factor',1.),
                "NumPart_ThisFile": self.num_part_this_file,
                "BoxSize": getattr(self,'box_size',0.),
                "Redshift": redshift,
                "Omega0": getattr(self,'omega_0',0.3),
                "OmegaLambda": getattr(self,'omega_lambda',0.7),
                "HubbleParam": getattr(self,'hubble_param',0.7),
                "Flag_Cooling": getattr(self,'flag_cooling',0),
                "Flag_StellarAge": getattr(self,'flag_stellar_age',0),
                "Flag_Sfr": getattr(self,'flag_sfr',0),
                "Flag_Metals": getattr(self,'flag_metals',0),
                "Flag_Feedback": getattr(self,'flag_feedback',0),
                "Flag_DoublePrecision": getattr(self,'flag_double_precision',0),
            })
        else:  # GADGET2/3 fallback
            redshift = 0.0 if self.scale_factor <= 0 else 1. / self.scale_factor - 1.
            header.attrs.update({
                "Time": getattr(self,'scale_factor',1.),
                "Redshift": redshift,
                "BoxSize": getattr(self,'box_size',0.),
                "Omega0": getattr(self,'omega_0',0.3),
                "OmegaLambda": getattr(self,'omega_lambda',0.7),
                "HubbleParam": getattr(self,'hubble_param',0.7),
                "NumPart_ThisFile": self.num_part_this_file,
                "Flag_Cooling": getattr(self,'flag_cooling',0),
                "Flag_StellarAge": getattr(self,'flag_stellar_age',0),
                "Flag_Sfr": getattr(self,'flag_sfr',0),
                "Flag_Metals": getattr(self,'flag_metals',0),
                "Flag_Feedback": getattr(self,'flag_feedback',0),
                "Flag_DoublePrecision": getattr(self,'flag_double_precision',0),
            })

        # Optional kwargs
        selection = getattr(self,'selection',None)
        if selection is not None:
            header.attrs.update({
                "Halo Centre": self.halo_centre,
                "Halo Systemic Velocity": self.halo_systemtic_velocity,
                "Halo Extent": self.halo_extent,
            })

        run_label = getattr(self,'run_label',None)
        
        if run_label is not None:
            header.attrs.update({
                "RunLabel": self.run_label,
            })
                
        header.attrs.update({
            "Periodic": getattr(self,'periodic',0),
        })

        if self.output_convention.upper() == 'SWIFT':
            cosmo = f.create_group('Cosmology')
            cosmo.attrs.update({
                "Omega_cdm": getattr(self,'omega_dm',0.255),
                "Omega_b": getattr(self,'omega_b',0.045),
                "Omega_lambda": getattr(self,'omega_lambda',0.7),
                "h": getattr(self,'hubble_param',0.7),
            })
        elif self.output_convention.upper() in ['GADGET4','AREPO']:
            params = f.create_group('Parameters')
            params.attrs.update({
                "BoxSize": getattr(self,'box_size',0.),
                "Omega0": getattr(self,'omega_0',0.3),
                "OmegaLambda": getattr(self,'omega_lambda',0.7),
                "HubbleParam": getattr(self,'hubble_param',0.7),
            })

    def _write_particles(self, f, **kwargs):
        for i, npart in enumerate(self.num_part_total):
            if npart > 0:
                group = f.create_group(f"PartType{i}")
                if getattr(self,'pos',None) is not None:
                    self._write_positions(group, i, self.idx, self.idx_type)
                if getattr(self,'vel',None) is not None:
                    self._write_velocities(group, i, self.idx, self.idx_type)
                if getattr(self,'pids',None) is not None:
                    self._write_ids(group, i, self.idx, self.idx_type)
                if getattr(self,'mass',None) is not None:
                    self._write_masses(group, i, self.idx, self.idx_type, kwargs.get("name_of_mass_block", "Masses"))
                if getattr(self,'u',None) is not None and i==getattr(self,'gas_type',0):
                    self._write_internal_energies(group, i, self.idx[self.idx_type[self.gas_type]:self.idx_type[self.gas_type+1]], self.idx_type, kwargs.get("name_of_u_block", "InternalEnergy"))
                if getattr(self,'rho',None) is not None and i==getattr(self,'gas_type',0):
                    self._write_densities(group, i, self.idx[self.idx_type[self.gas_type]:self.idx_type[self.gas_type+1]], self.idx_type)
                if getattr(self,'hsml',None) is not None:
                    if i==getattr(self,'gas_type',0):
                        self._write_smoothing_lengths(group, i, self.idx, self.idx_type)
                    elif i==getattr(self,'star_type',4):
                        self._write_smoothing_lengths(group, i, self.idx, self.idx_type)
                if getattr(self,'gas_Z',None) is not None and i==getattr(self,'gas_type',0):
                    self._write_metallicities(group, i, self.idx[self.idx_type[self.gas_type]:self.idx_type[self.gas_type+1]], self.idx_type,metallicity_type='gas')
                if getattr(self,'stellar_Z',None) is not None and i==getattr(self,'star_type',4):
                    self._write_metallicities(group, i, self.idx[self.idx_type[self.star_type]:self.idx_type[self.star_type+1]], self.idx_type,metallicity_type='stellar')
                if getattr(self,'sfr',None) is not None and i==getattr(self,'star_type',4):
                    self._write_star_formation_rates(group, i, self.idx[self.idx_type[self.star_type]:self.idx_type[self.star_type+1]], self.idx_type)
                if getattr(self,'age',None) is not None and i==getattr(self,'star_type',4):
                    self._write_stellar_ages(group, i, self.idx[self.idx_type[self.star_type]:self.idx_type[self.star_type+1]], self.idx_type)
                if getattr(self,'initmass',None) is not None and i==getattr(self,'star_type',4):
                    self._write_stellar_init_mass(group, i, self.idx[self.idx_type[self.star_type]:self.idx_type[self.star_type+1]], self.idx_type)

    def _write_positions(self, group, i, idx, idx_type):
        data_pos = group.create_dataset(
            "Coordinates",
            data=self.pos[idx][idx_type[i]:idx_type[i+1]],
        )
        data_pos.attrs.update({
            "CGSConversionFactor": self.unit_length_in_cgs,
            "aexp-scale-exponent": 1,
            "h-scale-exponent": -1,
        })

    def _write_velocities(self, group, i, idx, idx_type):
        data_vel = group.create_dataset(
            "Velocities",
            data=self.vel[idx][idx_type[i]:idx_type[i+1]],
        )
        data_vel.attrs.update({
            "CGSConversionFactor": self.unit_velocity_in_cgs,
            "aexp-scale-exponent": 0.5,
            "h-scale-exponent": 0,
        })
        
    def _write_ids(self, group, i, idx, idx_type):
        data_pids = group.create_dataset(
            "ParticleIDs",
            data=self.pids[idx][idx_type[i]:idx_type[i+1]],
        )
        data_pids.attrs.update({
            "CGSConversionFactor": 1,
            "aexp-scale-exponent": 0,
            "h-scale-exponent": 0,
        })

    def _write_masses(self, group, i, idx, idx_type, name_of_mass_block):
        data_mass = group.create_dataset(
            name_of_mass_block,
            data=self.mass[idx][idx_type[i]:idx_type[i+1]],
        )
        data_mass.attrs.update({
            "CGSConversionFactor": self.unit_mass_in_cgs,
            "aexp-scale-exponent": 0,
            "h-scale-exponent": -1,
        })

    def _write_internal_energies(self, group, i, idx, idx_type, name_of_u_block):
        data_u = group.create_dataset(
            name_of_u_block,
            data=self.u[idx][idx_type[i]:idx_type[i+1]],
        )
        data_u.attrs.update({
            "CGSConversionFactor": self.unit_velocity_in_cgs**2,
            "aexp-scale-exponent": 0,
            "h-scale-exponent": 0,
        })

    def _write_smoothing_lengths(self, group, i, idx, idx_type):
        data_smoothinglength = group.create_dataset(
            'SmoothingLength',
            data=self.hsml[idx][idx_type[i]:idx_type[i+1]],
        )
        data_smoothinglength.attrs.update({
            "CGSConversionFactor": self.unit_length_in_cgs,
            "aexp-scale-exponent": 1,
            "h-scale-exponent": -1,
        })

    def _write_densities(self, group, i, idx, idx_type):
        data_rho = group.create_dataset(
            'Density',
            data=self.rho[idx][idx_type[i]:idx_type[i+1]],
        )
        data_rho.attrs.update({
            "CGSConversionFactor": self.unit_density_in_cgs,
            "aexp-scale-exponent": 3,
            "h-scale-exponent": 2,
        })

    def _write_metallicities(self, group, i, idx, idx_type, metallicity_type='gas'):
        """Write gas or stellar metallicities into the HDF5 group."""
    
        if metallicity_type == 'gas':
            data = self.gas_metallicity[idx][idx_type[i]:idx_type[i+1]]
        elif metallicity_type == 'stellar':
            data = self.stellar_metallicity[idx][idx_type[i]:idx_type[i+1]]
        else:
            raise ValueError(f"Unknown metallicity_type '{metallicity_type}'")

        data_metal = group.create_dataset("Metallicity", data=data)
        data_metal.attrs.update({
            "CGSConversionFactor": 1,
            "aexp-scale-exponent": 0,
            "h-scale-exponent": 0,
        })

    def _write_star_formation_rates(self, group, i, idx, idx_type):
        data_sfr = group.create_dataset(
            'StarFormationRate',
            data=self.gas_sfr[idx][idx_type[i]:idx_type[i+1]],
        )
        data_sfr.attrs.update({
            "CGSConversionFactor": self.unit_sfr_in_cgs,
            "aexp-scale-exponent": 0,
            "h-scale-exponent": 0,
        })

    def _write_stellar_ages(self, group, i, idx, idx_type):
        data_age = group.create_dataset(
            'StellarFormationTime',
            data=self.stellarage[idx][idx_type[i]:idx_type[i+1]],
        )
        data_age.attrs.update({
            "CGSConversionFactor": 1,
            "aexp-scale-exponent": 0,
            "h-scale-exponent": 0,
        })

    def _write_stellar_init_mass(self, group, i, idx, idx_type):
        data_stellarinitmass = group.create_dataset(
            'StellarInitMass',
            data=self.stellarinitmass[idx][idx_type[i]:idx_type[i+1]],
        )
        data_stellarinitmass.attrs.update({
            "CGSConversionFactor": self.unit_mass_in_cgs,
            "aexp-scale-exponent": 1,
            "h-scale-exponent": -1,
        })

