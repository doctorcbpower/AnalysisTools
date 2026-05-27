import numpy as np
import os

class read_binary:
    """
    Requires the following attributes to be set before read_binary_snapshot():
    snaproot
    snapfileformat
    ismultifile
    gas_type
    star_type
    pids_type
    extra_blocks
    """

    def __init__(self):
        pass
    
    def read_binary_snapshot(self):
        """Main function to read binary snapshot data."""
        # Read header from the first file to get global particle counts
        self._read_binary_header(self.snapfilename)       

        idx_with_mass, NumPart, NumPart_InMassBlock = self._analyze_mass_structure()
        self._allocate_binary_memory(NumPart)

        blocknames, extra_flags = self._setup_extra_blocks_binary(NumPart)
        istart, ifinish = self._calculate_binary_offsets()

        # Iterate over the pre-validated file list, no reconstruction needed
        for filename in self.snapfiles:
            istart = self._read_single_binary_file(
                filename, istart, ifinish, idx_with_mass, blocknames, extra_flags
            )

    def _read_binary_header(self, filename):
        """Read header information from binary snapshot file."""
        print('Reading data from %s' % filename)
        self.NumPartType = 6

        with open(filename, 'rb') as f:
            offset = 0
            # Handle SNAP2 format offset - need to skip
            # initial 4 byte buffer, 4 chars, 1 int, 4 chars, 4 byte buffer.
            # This corresponds to 20 bytes in total, assuming a 4-byte int.
            if self.snapfileformat == 'SNAP2':
                offset += 20

            # Read header fields sequentially
            offset += 4
            f.seek(offset, os.SEEK_SET)
            self.num_part_this_file = np.fromfile(f, dtype=np.int32, count=6)

            offset += 24
            f.seek(offset, os.SEEK_SET)
            self.mass_table = np.fromfile(f, dtype=np.float64, count=6)
            
            offset += 48
            f.seek(offset, os.SEEK_SET)
            self.scale_factor = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 8
            f.seek(offset, os.SEEK_SET)
            self.redshift = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 16
            f.seek(offset, os.SEEK_SET)
            self.num_part_total = np.fromfile(f, dtype=np.int32, count=6)
            
            offset += 28
            f.seek(offset, os.SEEK_SET)
            self.num_files = np.fromfile(f, dtype=np.int32, count=1)[0]
            # FIX: don't let the header overwrite num_files — load_snapshot
            # already counted the actual files on disk, which is authoritative.
            # Only use the header value as a sanity check.
            if self.num_files != len(self.snapfiles):
                import warnings
                warnings.warn(
                    f"Header reports {self.num_files} file(s) but "
                    f"{len(self.snapfiles)} were found on disk. "
                    f"Using the on-disk count."
                )
            self.num_files = len(self.snapfiles)   # trust the filesystem
            
            offset += 4
            f.seek(offset, os.SEEK_SET)
            self.box_size = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 8
            f.seek(offset, os.SEEK_SET)
            self.omega_0 = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 8
            f.seek(offset, os.SEEK_SET)
            self.omega_lambda = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 8
            f.seek(offset, os.SEEK_SET)
            self.hubble_param = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            if self.num_files > 1:
                print('Data is split across %d files' % self.num_files)

    def _read_binary_header(self, filename):
        """Read header information from binary snapshot file."""
        print('Reading data from %s' % filename)

        self.NumPartType = 6  # GADGET standard
        
        with open(filename, 'rb') as f:
            offset = 0
  
            # Handle SNAP2 format offset - need to skip
            # initial 4 byte buffer, 4 chars, 1 int, 4 chars, 4 byte buffer.
            # This corresponds to 20 bytes in total, assuming a 4-byte int.
            if self.snapfileformat == 'SNAP2':
                offset += 20

            # Read header fields sequentially
            offset += 4
            f.seek(offset, os.SEEK_SET)
            self.num_part_this_file = np.fromfile(f, dtype=np.int32, count=6)

            offset += 24
            f.seek(offset, os.SEEK_SET)
            self.mass_table = np.fromfile(f, dtype=np.float64, count=6)
            
            offset += 48
            f.seek(offset, os.SEEK_SET)
            self.scale_factor = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 8
            f.seek(offset, os.SEEK_SET)
            self.redshift = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 16
            f.seek(offset, os.SEEK_SET)
            self.num_part_total = np.fromfile(f, dtype=np.int32, count=6)
            
            offset += 28
            f.seek(offset, os.SEEK_SET)
            self.num_files = np.fromfile(f, dtype=np.int32, count=1)[0]
            
            offset += 4
            f.seek(offset, os.SEEK_SET)
            self.box_size = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 8
            f.seek(offset, os.SEEK_SET)
            self.omega_0 = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 8
            f.seek(offset, os.SEEK_SET)
            self.omega_lambda = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            offset += 8
            f.seek(offset, os.SEEK_SET)
            self.hubble_param = np.fromfile(f, dtype=np.float64, count=1)[0]
            
            if self.num_files > 1:
                print('Data is split across %d files' % self.num_files)
    
    def _analyze_mass_structure(self):
        """Analyze which particle types are in the mass block."""
        idx_with_mass = np.where(self.mass_table == 0)[0]  # Species in mass block have zero MassTable entries
        
        NumPart = np.sum(self.num_part_total)
        print('Number of particles: %010d' % NumPart)
        
        NumPart_InMassBlock = np.sum(self.num_part_total[idx_with_mass])
        print('Number of particles in mass block: %010d' % NumPart_InMassBlock)
        
        return idx_with_mass, NumPart, NumPart_InMassBlock
    
    def _allocate_binary_memory(self, NumPart):
        """Allocate memory for particle data arrays."""
        self.pos = np.ndarray(shape=(NumPart, 3))
        self.vel = np.ndarray(shape=(NumPart, 3))
        
        if self.pids_type == 32:
            self.pids = np.ndarray(shape=(NumPart), dtype=np.uint32)
        else:
            self.pids = np.ndarray(shape=(NumPart), dtype=np.uint64)
        
        self.mass = np.ndarray(shape=(NumPart))
        
        # Gas particle specific arrays
        if self.num_part_total[0] > 0:
            self.u = np.ndarray(shape=(self.num_part_total[0]))
            self.rho = np.ndarray(shape=(self.num_part_total[0]))

    def _setup_extra_blocks_binary(self, NumPart):
        """Setup extra data blocks for binary format."""
        blocknames = self.extra_blocks
        extra_flags = {
            'isstellarage': False,
            'ismetallicity': False,
            'issfr': False,
            'ispotential': False
        }
        
        if len(self.extra_blocks) > 0:
            print('Loading extra blocks: %s' % self.extra_blocks)
            
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
            
            if 'POT' in self.extra_blocks:
                self.potential = np.ndarray(shape=(NumPart))
                extra_flags['ispotential'] = True
        
        return blocknames, extra_flags
    
    def _calculate_binary_offsets(self):
        """Calculate particle type offsets for binary format."""
        istart = np.zeros(self.NumPartType, dtype=np.uint64)
        istart[1:] = np.cumsum(self.num_part_total[:-1])
        ifinish = np.copy(istart)
        
        return istart, ifinish
    
    def _get_block_offset(self, base_offset=264):
        """Get initial block offset accounting for format differences."""
        offset = base_offset  # includes 2 x 4 byte buffers
        
        if self.snapfileformat == 'SNAP2':
            offset += 16
        
        return offset
    
    def _read_position_block(self, f, offset, NumPartInFile):
        """Read position data block."""
        offset += 4  # 1st 4 byte buffer
        
        if self.snapfileformat == 'SNAP2':
            offset += 16
        
        f.seek(offset, os.SEEK_SET)
        pos_block = np.fromfile(f, dtype=np.float32, count=3 * NumPartInFile)
        
        # Increment beyond the POS block
        offset += 3 * NumPartInFile * 4
        offset += 4  # 2nd 4 byte buffer
        
        return pos_block, offset
    
    def _read_velocity_block(self, f, offset, NumPartInFile):
        """Read velocity data block."""
        offset += 4  # 1st 4 byte buffer
        if self.snapfileformat == 'SNAP2':
            offset += 16
        
        f.seek(offset, os.SEEK_SET)
        vel_block = np.fromfile(f, dtype=np.float32, count=3 * NumPartInFile)
        
        # Increment beyond the VEL block
        offset += 3 * NumPartInFile * 4
        offset += 4  # 2nd 4 byte buffer
        
        return vel_block, offset
    
    def _read_ids_block(self, f, offset, NumPartInFile):
        """Read particle IDs block."""
        offset += 4  # 1st 4 byte buffer
        if self.snapfileformat == 'SNAP2':
            offset += 16
        
        f.seek(offset, os.SEEK_SET)
        
        if self.pids_type == 32:
            pids_block = np.fromfile(f, dtype=np.uint32, count=NumPartInFile)
            num_bytes = 4
        else:
            pids_block = np.fromfile(f, dtype=np.uint64, count=NumPartInFile)
            num_bytes = 8
        
        # Increment beyond the IDs block
        offset += NumPartInFile * num_bytes
        offset += 4  # 2nd 4 byte buffer
        
        return pids_block, offset
    
    def _read_mass_block(self, f, offset, NumPart_InMassBlock_InFile):
        """Read mass data block."""
        mass_block = None
        
        if NumPart_InMassBlock_InFile > 0:
            offset += 4  # 1st 4 byte buffer
            if self.snapfileformat == 'SNAP2':
                offset += 16
            
            f.seek(offset, os.SEEK_SET)
            mass_block = np.fromfile(f, dtype=np.float32, count=NumPart_InMassBlock_InFile)
            
            # Increment beyond the mass block
            offset += NumPart_InMassBlock_InFile * 4
            offset += 4  # 2nd 4 byte buffer
        
        return mass_block, offset
    
    def _read_gas_blocks(self, f, offset, NumPartInThisFile):
        """Read gas-specific data blocks (internal energy and density)."""
        u_block = None
        rho_block = None
        
        if NumPartInThisFile[0] > 0:
            # Read internal energy block
            offset += 4  # 1st 4 byte buffer
            if self.snapfileformat == 'SNAP2':
                offset += 16
            
            f.seek(offset, os.SEEK_SET)
            u_block = np.fromfile(f, dtype=np.float32, count=NumPartInThisFile[self.gas_type])
            
            # Increment beyond the internal energy block
            offset += NumPartInThisFile[0] * 4
            offset += 4  # 2nd 4 byte buffer
            
            # Read density block
            offset += 4  # 1st 4 byte buffer
            if self.snapfileformat == 'SNAP2':
                offset += 16
            
            f.seek(offset, os.SEEK_SET)
            rho_block = np.fromfile(f, dtype=np.float32, count=NumPartInThisFile[self.gas_type])
        
        return u_block, rho_block, offset
    
    def _read_extra_data_blocks(self, f, blocknames, extra_flags, NumPartInThisFile):
        """Read extra data blocks (metallicity, age, SFR, potential)."""
        extra_data = {}
        
        if extra_flags['ismetallicity']:
            offset = blocknames['Z'] + 20
            f.seek(offset, os.SEEK_SET)
            extra_data['gas_metals'] = np.fromfile(f, dtype=np.float32, count=NumPartInThisFile[self.gas_type])
            
            offset += 4 * NumPartInThisFile[0]
            f.seek(offset, os.SEEK_SET)
            extra_data['stellar_metals'] = np.fromfile(f, dtype=np.float32, count=NumPartInThisFile[self.star_type])
        
        if extra_flags['isstellarage']:
            offset = blocknames['AGE'] + 20
            f.seek(offset, os.SEEK_SET)
            extra_data['stellarage'] = np.fromfile(f, dtype=np.float32, count=NumPartInThisFile[self.star_type])
        
        if extra_flags['issfr']:
            offset = blocknames['SFR'] + 20
            f.seek(offset, os.SEEK_SET)
            extra_data['gas_sfr'] = np.fromfile(f, dtype=np.float32, count=NumPartInThisFile[self.gas_type])
        
        if extra_flags['ispotential']:
            offset = blocknames['POT'] + 20
            f.seek(offset, os.SEEK_SET)
            NumPartInFile = np.sum(NumPartInThisFile)
            extra_data['potential'] = np.fromfile(f, dtype=np.float32, count=NumPartInFile)
        
        return extra_data
    
    def _distribute_particle_data(self, istart, ifinish, NumPartInThisFile, idx_with_mass,
                                pos_block, vel_block, pids_block, mass_block, u_block, rho_block,
                                extra_data, extra_flags):
        """Distribute read data into particle type arrays."""
        astart = 0  # Position/velocity block index
        bstart = 0  # Single particle block index
        cstart = 0  # Mass block index
        
        for itype in range(self.NumPartType):
            if NumPartInThisFile[itype] > 0:
                afinish = astart + 3 * NumPartInThisFile[itype]
                bfinish = bstart + NumPartInThisFile[itype]
                
                # Positions and velocities
                self.pos[istart[itype]:ifinish[itype]] = pos_block[astart:afinish].reshape(NumPartInThisFile[itype], 3)
                self.vel[istart[itype]:ifinish[itype]] = vel_block[astart:afinish].reshape(NumPartInThisFile[itype], 3)
                self.pids[istart[itype]:ifinish[itype]] = pids_block[bstart:bfinish]
                
                # Masses
                if np.isin(itype, idx_with_mass):
                    cfinish = cstart + NumPartInThisFile[itype]
                    self.mass[istart[itype]:ifinish[itype]] = mass_block[cstart:cfinish]
                    cstart = cfinish
                else:
                    self.mass[istart[itype]:ifinish[itype]] = self.mass_table[itype] * np.ones(NumPartInThisFile[itype])
                
                # Potential
                if extra_flags['ispotential']:
                    self.potential[istart[itype]:ifinish[itype]] = extra_data['potential'][bstart:bfinish]
                
                # Gas data
                if NumPartInThisFile[0] > 0 and itype == self.gas_type:
                    self.u[istart[itype]:ifinish[itype]] = u_block[bstart:bfinish]
                    if rho_block.size > 0:
                        self.rho[istart[itype]:ifinish[itype]] = rho_block[bstart:bfinish]
                    
                    if extra_flags['ismetallicity']:
                        self.gas_metallicity[istart[itype]:ifinish[itype]] = extra_data['gas_metals'][bstart:bfinish]
                    
                    if extra_flags['issfr']:
                        self.gas_sfr[istart[itype]:ifinish[itype]] = extra_data['gas_sfr'][bstart:bfinish]
                
                # Star data
                if NumPartInThisFile[self.star_type] > 0 and itype == self.star_type:
                    if extra_flags['ismetallicity']:
                        self.stellar_metallicity[istart[itype]:ifinish[itype]] = extra_data['stellar_metals'][bstart:bfinish]
                    
                    if extra_flags['isstellarage']:
                        self.stellarage[istart[itype]:ifinish[itype]] = extra_data['stellarage'][bstart:bfinish]
                
                astart = afinish
                bstart = bfinish
    
    def _read_single_binary_file(self, filename, istart, ifinish, idx_with_mass, blocknames, extra_flags):
        """Read particle data from a single binary file."""
        with open(filename, 'rb') as f:
            # Read header to get particle counts for this file
            offset = 0
            if self.snapfileformat == 'SNAP2':
                offset += 16
            offset += 4
            f.seek(offset, os.SEEK_SET)
            NumPartInThisFile = np.fromfile(f, dtype=np.int32, count=6)
            
            NumPartInFile = np.sum(NumPartInThisFile)
            NumPart_InMassBlock_InFile = np.sum(NumPartInThisFile[idx_with_mass])
            
            if self.num_files > 1:
                print('Reading %010d particles from %s' % (NumPartInFile, filename))
                print('Number of particles in mass block: %010d' % NumPart_InMassBlock_InFile)
            
            # Start reading data blocks
            offset = self._get_block_offset()
            
            # Read main data blocks
            pos_block, offset = self._read_position_block(f, offset, NumPartInFile)
            vel_block, offset = self._read_velocity_block(f, offset, NumPartInFile)
            pids_block, offset = self._read_ids_block(f, offset, NumPartInFile)
            mass_block, offset = self._read_mass_block(f, offset, NumPart_InMassBlock_InFile)
            u_block, rho_block, offset = self._read_gas_blocks(f, offset, NumPartInThisFile)
            
            # Read extra data blocks
            extra_data = self._read_extra_data_blocks(f, blocknames, extra_flags, NumPartInThisFile)
            
            # Update finish indices for this file
            ifinish = istart + NumPartInThisFile.astype(np.uint64)
            
            # Distribute data to arrays
            self._distribute_particle_data(istart, ifinish, NumPartInThisFile, idx_with_mass,
                                         pos_block, vel_block, pids_block, mass_block,
                                         u_block, rho_block, extra_data, extra_flags)
            
            return ifinish
    
 
