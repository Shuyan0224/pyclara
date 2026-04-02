import numpy as np
import h5py
from scipy.constants import c, e, m_e, m_p
from fbpic.main import Simulation
from fbpic.openpmd_diag import FieldDiagnostic, ParticleDiagnostic
from fbpic.lpa_utils.bunch import add_particle_bunch_gaussian
from fbpic.lpa_utils.bunch import add_particle_bunch_from_arrays

class Fbpic_runner:
    def __init__(self):
        # --------------------
        # Simulation parameters
        # --------------------
        self.use_cuda = True
        self.n_order = -1  # Infinite-order solver
        self.Nm = 2        # Azimuthal modes
        self.write_dir = None

        # Plasma parameters
        self.n0 = 1e22  # Background plasma density [m^-3] (example value)
        # Calculate plasma wavelength
        self.omega_p = np.sqrt(self.n0 * e**2 / (m_e * 8.854e-12))
        self.lambda_p = 2 * np.pi * c / self.omega_p

        # Moving window domain size (meters)
        # Window size: 3*lambda_p long, 1*lambda_p radius
        self.zmax = 1.5 * self.lambda_p  # Moving window length
        self.zmin = -1.5 * self.lambda_p  # Start window centered around z=0
        self.rmax = 1 * self.lambda_p  # Radial extent = 1 plasma wavelength
        self.Nz = 512
        self.Nr = 96
        self.dz = (self.zmax - self.zmin) / self.Nz
        self.dr = self.rmax / self.Nr

        # More conservative time step
        self.dt = 0.5 * self.dz / c  # CFL condition for electromagnetic codes

        # Simulation time
        # Total propagation distance: 20 cm
        self.total_propagation = 2e-3  # 20 cm
        self.t_final = self.total_propagation / c  # Time for window to travel 20 cm
        self.N_steps = int(self.t_final / self.dt)

        # --------------------
        # Plasma parameters
        # --------------------
        # Plasma extends throughout the 20 cm simulation region
        # It will be continuously loaded as the window moves
        self.p_zmin = 1 * self.lambda_p
        self.p_zmax = self.total_propagation  # Fill the moving window
        self.p_rmax = self.rmax  # Fill radial extent (1*lambda_p)
        self.p_nz = 2   # Further reduced particles
        self.p_nr = 2
        self.p_nt = 4

        # --------------------
        # Beam parameters
        # --------------------
        # More reasonable beam energy (adjust as needed)
        self.total_charge = 250  # in pC
        self.gamma_b = None  
        self.beam_sigma_z = None  # Longitudinal RMS size 
        self.beam_sigma_r = None  # Transverse RMS size 
        self.beam_z0 = None    # Start before plasma

        # --------------------
        # Moving window and diagnostic 
        # --------------------
        self.v_window = c
        self.diag_period = 500

        # --------------------
        # Density profiles
        # --------------------
        def plasma_dens(z, r):
            """
            Plasma density profile with r^2 weighting for uniform 3D particle distribution

            In cylindrical coordinates, dV = r·dr·dθ·dz, so at larger r, volume increases.
            To compensate and get uniform particle density in 3D space, we scale n ∝ r^2.
            This ensures particles don't become sparse at large r.
            """
            # Create boolean masks for the plasma region
            in_z = (z >= self.p_zmin) & (z <= self.p_zmax)
            in_r = (r <= self.p_rmax)

            # Combine masks - plasma exists where both conditions are true
            in_plasma = in_z & in_r

            # Base density with r^2 weighting
            # Avoid division by zero: use (r/p_rmax)^2 normalized scaling
            r_normalized = r / self.p_rmax
            n = np.where(in_plasma, self.n0 , 0.0)

            return n

        # --------------------
        # Initialize simulation
        # --------------------
        print("Initializing simulation(1)...")
        self.sim = Simulation(self.Nz, self.zmax, self.Nr, self.rmax, self.Nm, self.dt, zmin=self.zmin,
                                n_order= self.n_order, use_cuda=self.use_cuda,
                                boundaries={'z':'open', 'r':'reflective'},)

        self.electrons = self.sim.add_new_species(
            q=-e, m=m_e, n=self.n0,
            dens_func= None, 
            p_zmin=self.p_zmin,  # Start where plasma actually exists
            p_zmax=self.p_zmax, 
            p_rmax=self.p_rmax,
            p_nz=self.p_nz, p_nr=self.p_nr, p_nt=self.p_nt
        )

        # Plasma ions
        self.protons = self.sim.add_new_species(
            q=e, m=m_p, n=self.n0,
            dens_func= None, 
            p_zmin=self.p_zmin,  # Match electrons
            p_zmax=self.p_zmax, 
            p_rmax=self.p_rmax,
            p_nz=self.p_nz, p_nr=self.p_nr, p_nt=self.p_nt
        )

        self.beam = None
        self.x_beam = None
        self.y_beam = None
        self.z_beam = None
        self.px_beam = None
        self.py_beam = None
        self.pz_beam = None
        self.x_beam = None
        self.n_particles = None
        self.beam_z0_off = None
        self.w_beam = None
        self.sim.set_moving_window(v=self.v_window)


    def set_input_Gaussian(self):
        self.gamma_b = self.total_charge / 0.511  # Mev/0.511 Reduced from extremely high value
        # Make beam much more compact: ~10 µm instead of ~84 µm
        self.beam_sigma_z = 1.4e-5  # Longitudinal RMS size ~8 µm
        self.beam_sigma_r = 1.85e-5  # Transverse RMS size ~8 µm
        self.beam_z0 = 1/2 * self.lambda_p    # Start at 105.6 µm (before plasma)

        self.beam = add_particle_bunch_gaussian(
                self.sim, q=-e, m=m_e, gamma0=self.gamma_b,
                n_emit=4.426719e-6,  # Zero emittance (idealized beam)
                sig_r=self.beam_sigma_r, sig_z=self.beam_sigma_z, sig_gamma=6.05,
                n_physical_particles= 1.5604e9,
                n_macroparticles= 262144,
                tf=0,  # Injection time
                zf=self.beam_z0,  # Beam center position
                z_injection_plane= 1 * self.lambda_p,
                boost=None,)
    
    def set_input_dict(self, input):
        # Read particle data
        # Read positions and monmentum
        for k, v in input.items():
            if k == "x":
                self.x_beam = np.array(v)
            if k == "y":
                self.y_beam = np.array(v)
            if k == "z":
                self.z_beam = np.array(v)
            if k == "px":
                self.px_beam = np.array(v)
            if k == "py":
                self.py_beam = np.array(v)
            if k == "p":
                self.pz_beam = np.array(v)
        self.n_particles = len(self.x_beam)
        print(f"The max logitonial velocity {np.max(self.pz_beam)}")
        print(f"  Loaded {self.n_particles} particles from file")

        # Apply z-offset to beam positions
        self.beam_z0_off = np.mean(self.z_beam) - 0.5*self.lambda_p
        self.z_beam = self.z_beam - self.beam_z0_off
        print(f"  Applied z-offset: -{self.beam_z0_off} m")
        print(f"  New z_beam range: [{self.z_beam.min():.6f}, {self.z_beam.max():.6f}] m")

        # Add beam species to simulation
        # w must be an array of weights, one per particle
        self.w_beam = np.full(self.n_particles, int(self.total_charge*1e-12/(self.n_particles*e)))  # Each macroparticle has weight 5950
        print(f"weigth = {self.w_beam}")
        self.beam = add_particle_bunch_from_arrays(self.sim, q=-e, m=m_e, x=self.x_beam,
                                                    y=self.y_beam, z=self.z_beam, ux=self.px_beam, 
                                                    uy =self.py_beam, uz=self.pz_beam, w=self.w_beam,
                                                    z_injection_plane = 1 * self.lambda_p)
        

    def run(self, outputdir):
        self.write_dir = outputdir

        # --------------------
        # Diagnostics
        # -------------------

        self.sim.diags = [
            # Save only E-field components (not B-field)
            FieldDiagnostic(
                self.diag_period, self.sim.fld, comm=self.sim.comm,
                fieldtypes=["E", "rho"],  # Only E-field and charge density
                write_dir=self.write_dir
            ),
            
            # Save all particle species
            ParticleDiagnostic(self.diag_period, {
                "beam": self.beam,'electrons': self.electrons
                }, comm=self.sim.comm,
                write_dir=self.write_dir
            )
        ]
        # --------------------
        # Run simulation
        # --------------------
        print("Starting simulation...")
        try:
            self.sim.step(self.N_steps)
            print("Simulation completed successfully!")
        except Exception as e:
            print(f"Simulation failed with error: {e}")
            print("Try reducing the time step or beam density further.")
                
    def get_output(self):
        pass


