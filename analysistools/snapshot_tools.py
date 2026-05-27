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
import logging
from typing import Union, List, Optional, Dict, Any

from .snapio_hdf5 import read_hdf5, write_hdf5
from .snapio_binary import read_binary

class SnapshotTools:
    """
    Handle cosmological simulation snapshots (HDF5 or GADGET binary).
    
    Usage:
    # Configure once, then read later:
    snap = SnapshotTools(snapfileformat="HDF5")
    snap.load_snapshot("snap_010.hdf5")
    data = snap.read()
    
    # Or do it in one go:
    snap = SnapshotTools(snapfileformat="HDF5")
    data = snap.read("snap_010.hdf5")
    """
    
    def __init__(
        self,
        snapfileformat: Union[int, str] = "HDF5",
        gas_type: int = 0,
        dm_type: int = 1,
        star_type: int = 4,
        bh_type: int = 5,
        **kwargs,
    ):
        self.snapfilename: Optional[str] = None
        
        # --- format normalisation ---
        fmtmap = {
            1: "SNAP1", 2: "SNAP2", 3: "HDF5",
            "1": "SNAP1", "2": "SNAP2", "3": "HDF5"
        }
        self.snapfileformat = fmtmap.get(snapfileformat, str(snapfileformat).upper())
        if self.snapfileformat not in {"SNAP1", "SNAP2", "HDF5"}:
            raise ValueError(f"Unknown snapshot format: {snapfileformat}")
        
        # --- particle types ---
        self.gas_type = gas_type
        self.dm_type = dm_type
        self.star_type = star_type
        self.bh_type = bh_type
        
        # --- options ---
        self.convention = kwargs.get("convention", "GADGET2/3")
        self.positions_only = kwargs.get("positions_only", False)
        self.hires_only = kwargs.get("hires_only", False)
        self.get_ptypes = kwargs.get("get_ptypes", False)
        self.extra_blocks: List[str] = kwargs.get("extra_blocks", [])
        self.positions_type = kwargs.get("positions_type", "float32")
        self.pids_type = kwargs.get("pids_type", 32)
        self.not_hires_ptypes = kwargs.get("not_hires_ptypes", [2, 3, 7])
        self.isics = kwargs.get("isics", False)
        self.is_multifile = False
        self.binary_suffix = kwargs.get("binary_suffix",".dat")   # NEW
                
        # Set num_part_type for consistency
        self.num_part_type = 6
        
        self._set_units()

        # Logging setup
        self.logger = logging.getLogger(__name__)
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
            self.logger.addHandler(handler)
            self.logger.propagate = False
        self.logger.setLevel(kwargs.get("loglevel", logging.INFO))
        
    def _set_units(self):
        """Define useful units needed for snapshots"""
        SOLAR_MASS_IN_CGS = 1.989e33
        YEAR_IN_CGS = 60*60*24*365
        KPC_IN_CGS = 3.0856e21
        KM_PER_SEC_IN_CGS = 1.e5
        RHOCRIT0 = 27.755
        
        self.unit_length_in_cgs = KPC_IN_CGS  # kpc
        self.unit_mass_in_cgs = 1e10 * SOLAR_MASS_IN_CGS  # 1e10 Msol
        self.unit_velocity_in_cgs = KM_PER_SEC_IN_CGS  # km/s
        self.unit_time_in_cgs = self.unit_length_in_cgs / self.unit_velocity_in_cgs
        self.unit_density_in_cgs = self.unit_mass_in_cgs / self.unit_length_in_cgs**3
        self.unit_sfr_in_cgs = SOLAR_MASS_IN_CGS / YEAR_IN_CGS
    
    def load_snapshot(self, snapfilename: str) -> None:
        """Register and validate the snapshot file for later use."""
        import os
        import glob

        # Strip known suffixes to get a clean base
        base = snapfilename
        for suffix in (".hdf5", ".0.hdf5", ".0", self.binary_suffix):
            if base.endswith(suffix):
                base = base[: -len(suffix)]
                break  # only strip one suffix

        self.snaproot = base

        # --- probe filesystem ---
        hdf5_single = base + ".hdf5"
        hdf5_multi  = sorted(glob.glob(base + ".[0-9]*.hdf5"))
        bin_single  = base + self.binary_suffix                              # NEW: single binary
        bin_multi   = sorted(glob.glob(base + ".[0-9]*"))
        bin_multi   = [f for f in bin_multi if not f.endswith(".hdf5")]

        if os.path.exists(hdf5_single):
            self.snapfileformat = "HDF5"
            self.is_multifile   = False
            self.snapfilename   = hdf5_single
            self.snapfiles      = [hdf5_single]
            self.num_files      = 1

        elif hdf5_multi:
            self.snapfileformat = "HDF5"
            self.is_multifile   = True
            self.snapfilename   = hdf5_multi[0]
            self.snapfiles      = hdf5_multi
            self.num_files      = len(hdf5_multi)

        elif bin_single and os.path.isfile(bin_single):
            self.is_multifile = False
            self.snapfilename = bin_single
            self.snapfiles    = [bin_single]
            self.num_files    = 1
            if self.snapfileformat not in {"SNAP1", "SNAP2"}:
                self.snapfileformat = "SNAP2"

        elif bin_multi:
            self.is_multifile = True
            self.snapfilename = bin_multi[0]
            self.snapfiles    = bin_multi
            self.num_files    = len(bin_multi)
            if self.snapfileformat not in {"SNAP1", "SNAP2"}:
                self.snapfileformat = "SNAP2" 
        else:
            raise FileNotFoundError(
                f"No snapshot found for '{snapfilename}'. "
                f"Tried: '{hdf5_single}', multi-file HDF5, "
                f"'{bin_single}', and numbered binary files."
            )

        self.snapbase = base
        self.logger.info(
            f"Snapshot detected: format={self.snapfileformat}, "
            f"multifile={self.is_multifile}, "
            f"files={self.num_files}, "
            f"base='{self.snapbase}'"
        )
    
    def _transfer_attributes_to_reader(self, reader):
        """Transfer snapshot metadata to the reader."""

        required = {
            "snaproot": getattr(self, "snaproot", None),
            "snapfilename":  getattr(self, "snapfilename", None),  
            "snapfileformat": self.snapfileformat,
            "ismultifile": self.is_multifile,
            "snapfiles": getattr(self, "snapfiles", None),
        }

        optional = {
            "dm_type": self.dm_type,
            "gas_type": self.gas_type,
            "star_type": self.star_type,
            "bh_type": self.bh_type,
            "positions_type": self.positions_type,
            "positions_only": self.positions_only,
            "pids_type": self.pids_type,
            "extra_blocks": self.extra_blocks,
            "hires_only": self.hires_only,
            "get_ptypes": self.get_ptypes,
            "not_hires_ptypes": self.not_hires_ptypes,
            "num_part_type": self.num_part_type,
            "isics": self.isics,
        }

        # Required: must exist
        for name, value in required.items():
            if value is None:
                raise AttributeError(
                    f"SnapshotTools missing required attribute '{name}' for reader"
                )
            setattr(reader, name, value)

        # Optional: override reader defaults
        for name, value in optional.items():
            if value is not None:
                setattr(reader, name, value)
    
    def _transfer_attributes_to_writer(self, writer):
        """Transfer configuration attributes to the writer instance."""
        attrs_to_transfer = [
            'gas_type', 'dm_type', 'star_type', 'bh_type', 'convention',
            'positions_only', 'hires_only', 'extra_blocks', 'positions_type',
            'pids_type', 'not_hires_ptypes', 'name_of_mass_block', 'name_of_u_block',
            'unit_length_in_cgs', 'unit_mass_in_cgs', 'unit_velocity_in_cgs',
            'unit_time_in_cgs', 'unit_density_in_cgs', 'unit_sfr_in_cgs',
            'num_part_type', 'num_part_total', 'num_part_this_file', 'mass_table',
            'box_size', 'dimension', 'scale_factor', 'redshift', 'time',
            'omega_dm', 'omega_b', 'omega_0', 'omega_lambda',
            'h', 'hubble_param',
            'flag_cooling', 'flag_stellar_age', 'flag_sfr', 'flag_metals',
            'flag_feedback', 'flag_double_precision',
            'selection', 'halo_centre', 'halo systemic_velocity', 'halo_extent',
            'run_label', 'periodic',
        ]
        
        for attr in attrs_to_transfer:
            if hasattr(self, attr):
                setattr(writer, attr, getattr(self, attr))

    def _transfer_datasets_to_writer(self, writer, datasets_to_transfer, copy=False):
        """Transfer datasets to the writer, with optional deep copy."""
        for dset in datasets_to_transfer:
            if hasattr(self, dset):
                data = getattr(self, dset)
                if copy:
                    if hasattr(data, 'shape'):  # ndarray or h5py.Dataset
                        data = np.array(data)   # make a real copy
                    else:
                        # fallback for generic mutable containers
                        import copy as pycopy
                        data = pycopy.deepcopy(data)
                setattr(writer, dset, data)
    def _initialize_writer_header(self, writer, idx, idx_type):
        # --- particle counts ---
        num_part_this_file = np.zeros(6, dtype=np.int64)

        for ptype in range(6):
            num_part_this_file[ptype] = np.sum(idx_type == ptype)

        writer.num_part_this_file = num_part_this_file

        # Often required as well
        writer.num_part_total = num_part_this_file.copy()

        # --- mass table (0 = individual masses stored) ---
        writer.mass_table = np.zeros(6, dtype=np.float64)

        # --- basic cosmology / box ---
        writer.scale_factor = getattr(self, "scale_factor", 1.0)
        writer.box_size = getattr(self, "box_size", 0.0)

        writer.omega_0 = getattr(self, "omega_0", 0.3)
        writer.omega_lambda = getattr(self, "omega_lambda", 0.7)
        writer.hubble_param = getattr(self, "hubble_param", 0.7)

        # --- flags ---
        writer.flag_double_precision = 1
        writer.flag_cooling = 0
        writer.flag_sfr = 0
        writer.flag_stellar_age = 0
        writer.flag_metals = 0
        writer.flag_feedback = 0

    def _transfer_attributes_from_reader(self, reader):
        """Transfer data attributes from reader back to main instance."""
        # Transfer all attributes that don't start with underscore
        for attr_name in dir(reader):
            if not attr_name.startswith('_') and not callable(getattr(reader, attr_name)):
                setattr(self, attr_name, getattr(reader, attr_name))
    
    def read_snapshot(self, filename: Optional[str] = None, convention: Optional[str] = None):
        """
        Read the snapshot data.
        
        Parameters
        ----------
        filename : str, optional
            If given, sets/overrides the snapshot file to read.
        convention : str, optional
            If given, overrides the simulation convention.
            
        Returns
        -------
        self : SnapshotTools
            Returns self for method chaining, with data loaded as attributes.
        """
        print("Reading snapshot data...")
       
        if filename is not None:
            self.load_snapshot(filename)
        
        if convention is not None:
            self.convention = convention

        try:
            self.logger.info(f"Reading snapshot '{filename}' using {self.snapfileformat}")
            if self.snapfileformat == "HDF5":
                reader = read_hdf5()
                self._transfer_attributes_to_reader(reader)
                reader.read_hdf5_snapshot(snapfilename=self.snapfilename,convention=self.convention)
                self._transfer_attributes_from_reader(reader)
            else:
                reader = read_binary()
                self._transfer_attributes_to_reader(reader)
                reader.read_binary_snapshot()
                self._transfer_attributes_from_reader(reader)
                
        except Exception as e:
            logging.error(f"Error reading snapshot {self.snapfilename}: {e}")
            raise
        
        return self
        
    def write_snapshot(self, filename: str,
                             idx: np.int32,
                             idx_type: np.int64,
                             file_format: str = "hdf5",
                             convention: str = "SWIFT",
                             blocks_to_write: list[str]| None = None,
                             **kwargs
                       ) -> None:
        """
        Write snapshot data to file in the chosen format.

        Parameters
        ----------
        filename : str
            Path to the output file.
        file_format : str, optional
            File format to use. Currently supports "hdf5".
        convention: str, optional
            Code convention to use - determines the names of blocks
        blocks_to_write : list[str], optional   
            List of dataset names to include in the output. Defaults to
            ['pos', 'vel', 'pids', 'mass', 'groupid'].
        Optional kwargs can include metadata like:
        halo_centre, halo_systemic_velocity, halo_extent, run_label, etc.
        """
        blocks = blocks_to_write or ['pos', 'vel', 'pids', 'mass', 'groupid']

        self.logger.info(f"Writing snapshot '{filename}' using {self.snapfileformat}")

        if file_format.lower() == "hdf5":
            writer = write_hdf5(output_convention=convention, idx=idx, idx_type=idx_type, **kwargs)

            self._transfer_attributes_to_writer(writer)

            self._initialize_writer_header(writer, idx=idx, idx_type=idx_type)
            """
            Options are: pos, vel, pids, mass, u, rho, hsml, gas_Z, stellar_Z,
                         sfr, age, initmass, groupid
            """

            self._transfer_datasets_to_writer(writer,
                                              blocks,
                                              )

            writer.write_hdf5_snapshot(filename)
        else:
            raise ValueError(f"Unsupported snapshot format: {file_format}")

    def ParticleOffsetsByType(self,NumPartByType):
        """
        Return offsets of particle species by type
        """
        return np.concatenate([[0],np.cumsum(NumPartByType)],dtype=np.int64)
    
    def LoadParticlesByType(self, part_type: str = 'all'):
        """
        Load particles into separate objects by type.
        
        Parameters
        ----------
        part_type : str
            Which particle types to load. Options: 'all', 'gas', 'star', 'dm', 'bh'
        """
        if not hasattr(self, 'pos') or not hasattr(self, 'num_part_total'):
            raise RuntimeError("No data loaded. Call read() first.")
        
        # Calculate offsets for each particle type
        offsets = self._calculate_particle_offsets()
        
        # Determine which types to load
        load_flags = self._determine_load_flags(part_type, offsets)
        
        # Initialize potential array if needed
        if not hasattr(self, 'potential') or self.potential is None:
            self.potential = np.zeros(shape=(np.sum(self.num_part_total)), dtype=np.float32)
        
        # Load particle data for each requested type
        if load_flags['gas']:
            print('gas',offsets['gas'])
            self.gas = self._create_particle_object('gas', offsets)
        
        if load_flags['dm']:
            print('dm',offsets['dm'])
            self.dm = self._create_particle_object('dm', offsets)
        
        if load_flags['star']:
            print('star',offsets['star'])
            self.star = self._create_particle_object('star', offsets)
        
        if load_flags['bh']:
            print('bh',offsets['bh'])
            self.bh = self._create_particle_object('bh', offsets)
    
    def _calculate_particle_offsets(self):
        """Calculate start/end offsets for each particle type."""
        offsets = {}
        
        offsets['gas'] = {
            'start': np.sum(self.num_part_total[:self.gas_type]),
            'end': np.sum(self.num_part_total[:self.gas_type + 1])
        }
        offsets['dm'] = {
            'start': np.sum(self.num_part_total[:self.dm_type]),
            'end': np.sum(self.num_part_total[:self.dm_type + 1])
        }
        offsets['star'] = {
            'start': np.sum(self.num_part_total[:self.star_type]),
            'end': np.sum(self.num_part_total[:self.star_type + 1])
        }
        offsets['bh'] = {
            'start': np.sum(self.num_part_total[:self.bh_type]),
            'end': np.sum(self.num_part_total[:self.bh_type + 1])
        }
        
        return offsets
    
    def _determine_load_flags(self, part_type: str, offsets: Dict):
        """Determine which particle types should be loaded."""
        load_flags = {'gas': False, 'dm': False, 'star': False, 'bh': False}
        
        if part_type == 'all':
            for ptype in ['gas', 'star', 'bh', 'dm']:
                if offsets[ptype]['end'] - offsets[ptype]['start'] > 0:
                    load_flags[ptype] = True
        elif part_type in load_flags:
            if offsets[part_type]['end'] - offsets[part_type]['start'] > 0:
                load_flags[part_type] = True
        else:
            raise ValueError(f"Unknown particle type: {part_type}")
        
        return load_flags
    
    def _create_particle_object(self, ptype: str, offsets: Dict):
        """Create a ParticleProperties object for the specified type."""
        start, end = offsets[ptype]['start'], offsets[ptype]['end']
        
        kwargs = {}
        if ptype == 'gas':
            if hasattr(self, 'u'):
                kwargs['internal_energy'] = self.u[start:end]
            if hasattr(self, 'rho'):
                kwargs['density'] = self.rho[start:end]
        
        if hasattr(self, 'groupid'):
            kwargs['groupid'] = self.groupid[start:end]
        
        return self.ParticleProperties(
            self.pos[start:end],
            self.vel[start:end],
            self.pids[start:end],
            self.mass[start:end],
            self.potential[start:end],
            **kwargs
        )
    
    class ParticleProperties:
        """Container for particle data of a specific type."""
        
        def __init__(self, pos, vel, pids, mass, potential, **kwargs):
            self.pos = pos
            self.vel = vel
            self.pids = pids
            self.mass = mass
            self.potential = potential
            
            # Optional gas-specific properties
            if 'internal_energy' in kwargs:
                self.internal_energy = kwargs['internal_energy']
            if 'density' in kwargs:
                self.density = kwargs['density']
            if 'groupid' in kwargs:
                self.groupid = kwargs['groupid']
    
    def UnitConversion(self, **kwargs):
        """
        Apply unit conversions to the data.
        
        Parameters
        ----------
        convert_to_physical : bool, optional
            Convert from comoving to physical coordinates
        convert_to_comoving : bool, optional
            Convert from physical to comoving coordinates
        convert_to_per_littleh : bool, optional
            Convert to units per little h
        convert_to_littleh : bool, optional
            Convert from units per little h
        """
        if not hasattr(self, 'pos') or not hasattr(self, 'scale_factor'):
            raise RuntimeError("No data loaded or missing cosmological parameters")
        
        if kwargs.get('convert_to_physical'):
            self.pos *= self.ScaleFactor
            if hasattr(self, 'BoxSize'):
                self.BoxSize *= self.ScaleFactor
        
        if kwargs.get('convert_to_comoving'):
            self.pos /= self.ScaleFactor
            if hasattr(self, 'BoxSize'):
                self.BoxSize /= self.ScaleFactor
        
        if kwargs.get('convert_to_per_littleh'):
            self.pos *= self.HubbleParam
            if hasattr(self, 'mass'):
                self.mass *= self.HubbleParam
            if hasattr(self, 'BoxSize'):
                self.BoxSize *= self.HubbleParam
        
        if kwargs.get('convert_to_littleh'):
            self.pos /= self.HubbleParam
            if hasattr(self, 'mass'):
                self.mass /= self.HubbleParam
            if hasattr(self, 'BoxSize'):
                self.BoxSize /= self.HubbleParam


# Utility functions
def select_particles(val, valoffset, size, geometry, **kwargs):
    """
    Select particles within a specified region.
    
    Parameters
    ----------
    val : array_like
        Particle positions
    valoffset : array_like
        Center of selection region
    size : float
        Size of selection region (box size or radius)
    geometry : str
        Selection geometry ('cubic' or 'spherical')
    **kwargs : dict
        Additional options (periodic, scale_length, etc.)
    
    Returns
    -------
    ipick : array_like
        Indices of selected particles
    """
    dval = val - valoffset
    
    # Handle periodicity
    if kwargs.get('periodic') and kwargs.get('scale_length') is not None:
        scale_length = kwargs['scale_length']
        dval = np.where(dval > 0.5 * scale_length, dval - scale_length, dval)
        dval = np.where(dval < -0.5 * scale_length, dval + scale_length, dval)
    else:
        logging.info('Ignoring periodicity')
    
    # Apply geometric selection
    if geometry == 'cubic':
        ipick = np.logical_and(np.abs(dval[:, 0]) < size/2, np.abs(dval[:, 1]) < size/2)
        ipick = np.logical_and(ipick, np.abs(dval[:, 2]) < size/2)
    elif geometry == 'spherical':
        r2 = dval[:, 0]**2 + dval[:, 1]**2 + dval[:, 2]**2
        ipick = np.where(r2 < size*size)[0]
    else:
        raise ValueError(f'Unknown geometry: {geometry}')
    
    if kwargs.get('return_ptype'):
        return ipick, kwargs.get('ptype')[ipick]
    else:
        return ipick


def place_points_in_mesh(pos, pos_offset, size, mesh_dimension, **kwargs):
    """
    Place particles into a mesh grid.
    
    Parameters
    ----------
    pos : array_like
        Particle positions
    pos_offset : array_like
        Offset for mesh origin
    size : float
        Size of mesh region
    mesh_dimension : int
        Number of mesh cells per dimension
    
    Returns
    -------
    array_like
        Mesh coordinates for each particle
    """
    return np.fix(mesh_dimension * (pos - pos_offset) / size)

 
