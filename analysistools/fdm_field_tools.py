"""
fdm_field_tools.py

Mesh-field analog to SnapshotTools, for reading FDM wavefunction field
snapshots (src/fdm/fdm_io.c's format in arepo_solas) -- NOT particle
data, so it does not fit into SnapshotTools/LoadParticlesByType at all.
Styled to match this package's existing construction pattern (a
lightweight constructor, then simple attribute access) rather than
introducing a different interface style for field data specifically.

File format read (see arepo_solas src/fdm/fdm_io.c for the authoritative
spec): HDF5, attributes "N" (int), "L" (double, box size), "FDMMass"
(double, mc^2 in eV); datasets "PsiReal"/"PsiImag", each N*N*N doubles
in plain global (x,y,z) row-major order (x slowest, z fastest).
"""

import numpy as np
import h5py


class FDMFieldTools:
    """
    Loads a single FDM wavefunction field snapshot.

    Usage mirrors SnapshotTools: construct, then either read immediately
    or defer:

        field = FDMFieldTools()
        field.load_field("fdm_field_010.hdf5")
        field.read_field()

    or combined:

        field = FDMFieldTools("fdm_field_010.hdf5")

    After reading, data is available as attributes:

        field.psi       # complex128 array, shape (N, N, N)
        field.density   # |psi|^2, real array, shape (N, N, N)
        field.N         # grid size
        field.L         # box size (code length units)
        field.mass      # mc^2, eV
    """

    def __init__(self, filename=None):
        self.filename = None
        self.N = None
        self.L = None
        self.mass = None
        self.psi = None

        if filename is not None:
            self.load_field(filename)
            self.read_field()

    def load_field(self, filename):
        """Registers the file to read -- mirrors SnapshotTools.load_snapshot's
        separation of "register" from "read", even though for a single
        field file (unlike multi-file particle snapshots) there is no
        file-list detection to do here."""
        self.filename = filename

    def read_field(self, filename=None):
        if filename is not None:
            self.filename = filename
        if self.filename is None:
            raise ValueError("No filename set -- call load_field() first or pass one to read_field().")

        with h5py.File(self.filename, "r") as f:
            self.N = int(f.attrs["N"])
            self.L = float(f.attrs["L"])
            self.mass = float(f.attrs["FDMMass"])
            psi_real = f["PsiReal"][:]
            psi_imag = f["PsiImag"][:]

        self.psi = psi_real + 1j * psi_imag
        return self

    @property
    def density(self):
        """|psi|^2 -- NOT yet multiplied by the boson mass; this is the
        quantity fdm_update_potential() itself calls "dens" before its
        own m_code multiplication (see arepo_solas src/fdm/fdm_poisson.c).
        Multiply by mass (code units) separately if a physical mass
        density is specifically needed, rather than baking a unit
        system assumption into this general-purpose reader."""
        if self.psi is None:
            raise RuntimeError("No field loaded -- call read_field() first.")
        return np.abs(self.psi) ** 2

    @property
    def dx(self):
        """Cell spacing, code length units."""
        return self.L / self.N

    def central_density_exact(self, center=None, quantity=None):
        """The density at the single grid cell nearest to `center`, NOT a
        coarse radial-bin average -- for a steeply-peaked profile (like a
        soliton), the first radial bin's average can substantially
        underestimate the true central value if bins are wide relative to
        the profile's own curvature near r=0, which then systematically
        biases anything measured relative to "central density" (e.g.
        half_density_radius in validate_soliton_stability.py) outward.
        Use this, not profile[0], as the reference value for any such
        measurement."""
        if quantity is None:
            quantity = self.density
        if center is None:
            center = np.array([self.L / 2.0] * 3)

        i = int(round(center[0] / self.dx))
        j = int(round(center[1] / self.dx))
        k = int(round(center[2] / self.dx))
        i, j, k = np.clip([i, j, k], 0, self.N - 1)
        return quantity[i, j, k]

    def radial_profile(self, quantity=None, center=None, nbins=50, rmax=None):
        """
        Spherically-averaged radial profile of `quantity` (defaults to
        self.density) around `center` (defaults to the box center,
        matching the convention used throughout arepo_solas's own FDM
        validation tests, e.g. fdm_test_poisson.c).

        Returns (r_bin_centers, profile_values, counts_per_bin) -- counts
        included so callers can judge shot-noise/binning reliability at
        small r, not just get a silently-noisy innermost bin.
        """
        if quantity is None:
            quantity = self.density
        if center is None:
            center = np.array([self.L / 2.0] * 3)
        if rmax is None:
            rmax = self.L / 2.0

        coords = (np.arange(self.N) + 0.0) * self.dx
        X, Y, Z = np.meshgrid(coords, coords, coords, indexing="ij")
        r = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2 + (Z - center[2]) ** 2)

        bin_edges = np.linspace(0, rmax, nbins + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        profile = np.zeros(nbins)
        counts = np.zeros(nbins, dtype=int)

        r_flat = r.ravel()
        q_flat = quantity.ravel()
        bin_idx = np.digitize(r_flat, bin_edges) - 1

        for b in range(nbins):
            mask = bin_idx == b
            counts[b] = np.sum(mask)
            if counts[b] > 0:
                profile[b] = np.mean(q_flat[mask])

        return bin_centers, profile, counts

    def total_norm(self):
        """Sum|psi|^2 over the whole field -- the conserved quantity
        checked throughout arepo_solas's own validation (fdm_test_step.c
        etc)."""
        return np.sum(self.density)

    def hbar_over_m_code(self, unit_length_cm=3.085678e21, unit_mass_g=1.989e43, unit_velocity_cm_s=1e5):
        """hbar/m in code units -- IDENTICAL formula to
        fdm_hbar_over_m_code() (src/fdm/fdm_integrator.c) and
        wavepacket_construction.py's own hbar_over_m_code() (used for
        the soliton+NFW halo IC's own wave-packet construction) --
        reproduced here a third time rather than imported, since this
        module is meant to be usable standalone within AnalysisTools
        without a dependency on the IC-generation code, but the
        CONSTANTS AND FORMULA must stay byte-for-byte identical across
        all three copies. Defaults match this project's own established
        kpc + 1e10 Msun + km/s unit convention -- override if a
        snapshot was generated with different units.
        """
        PLANCK = 6.6260695e-27  # erg*s (h, not hbar)
        CLIGHT = 2.99792458e10
        ELECTRONVOLT_IN_ERGS = 1.60217656e-12
        hbar = PLANCK / (2.0 * np.pi)
        m_grams = (self.mass * ELECTRONVOLT_IN_ERGS) / CLIGHT ** 2
        hbar_over_m_cgs = hbar / m_grams
        unit_time_s = unit_length_cm / unit_velocity_cm_s
        return hbar_over_m_cgs / (unit_length_cm ** 2 / unit_time_s)

    def density_slice(self, axis="z", index=None):
        """A 2D slice of self.density through the field at a fixed index
        along `axis` (defaults to the box-center index, matching this
        project's own box-centered halo/soliton convention). Returns
        (slice_2d, coord1, coord2) where coord1/coord2 are the 1D
        coordinate arrays for the two in-plane axes -- e.g. for
        axis="z", returns (density[:,:,index], x_coords, y_coords).
        Plot with e.g. imshow(np.log10(slice_2d.T), extent=[...]) for
        the standard density-slice visualization used throughout this
        project (log scale -- density spans many orders of magnitude
        between a soliton core and its surrounding envelope)."""
        if index is None:
            index = self.N // 2
        coords = (np.arange(self.N) + 0.5) * self.dx

        if axis == "z":
            return self.density[:, :, index], coords, coords
        elif axis == "y":
            return self.density[:, index, :], coords, coords
        elif axis == "x":
            return self.density[index, :, :], coords, coords
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got {axis!r}")

    def velocity_field_slice(self, axis="z", index=None, **hbar_m_kwargs):
        """The in-plane velocity field (probability current,
        v = (hbar/m)*Im(psi* grad(psi))/|psi|^2 -- see this project's
        own visualization work for the full reasoning) on a 2D slice
        through the field. Returns (vx, vy, coord1, coord2), all in
        km/s if the default unit convention applies (pass
        unit_velocity_cm_s etc via hbar_m_kwargs to override, matching
        hbar_over_m_code()'s own signature).

        Formally singular where density->0 (destructive-interference
        nodes/vortex lines) -- a genuine physical feature of the
        interference structure, not a numerical artifact; not clipped
        or smoothed away here, since doing so would hide exactly the
        structure this method exists to help visualize. Callers should
        expect occasional very large values near such nodes and handle
        them at the visualization stage (e.g. percentile-based color/
        arrow-length capping), not by treating them as errors.
        """
        if index is None:
            index = self.N // 2
        coords = (np.arange(self.N) + 0.5) * self.dx
        hbar_m = self.hbar_over_m_code(**hbar_m_kwargs)

        if axis == "z":
            psi_slice = self.psi[:, :, index]
        elif axis == "y":
            psi_slice = self.psi[:, index, :]
        elif axis == "x":
            psi_slice = self.psi[index, :, :]
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got {axis!r}")

        density_slice = np.abs(psi_slice) ** 2
        dpsi_d1 = np.gradient(psi_slice, self.dx, axis=0)
        dpsi_d2 = np.gradient(psi_slice, self.dx, axis=1)

        v1 = hbar_m * np.imag(np.conj(psi_slice) * dpsi_d1) / (density_slice + 1e-300)
        v2 = hbar_m * np.imag(np.conj(psi_slice) * dpsi_d2) / (density_slice + 1e-300)

        # hbar_over_m_code's own natural output is in code velocity
        # units (km/s under this project's default convention) already,
        # given hbar_m itself carries units of length^2/time and
        # dpsi_d/self.dx contributes an extra 1/length -- no further
        # conversion needed here.
        return v1, v2, coords, coords
