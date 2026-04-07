import numpy as np
import pyclara
import h5py
from scipy.constants import c, e, m_e, m_p, epsilon_0
from fbpic.main import Simulation
from fbpic.openpmd_diag import FieldDiagnostic, ParticleDiagnostic
from fbpic.lpa_utils.bunch import add_particle_bunch_gaussian
from fbpic.lpa_utils.bunch import add_particle_bunch_from_arrays



class Fbpic_runner:
    def __init__(self, n0 = 1e22):
        # Simulation parameters
        self.set_Sim_control()
        # Plasma parameters
        self.n0 = n0  # Background plasma density [m^-3] (example value)
        self.get_lambda_p()
        # Moving window domain size (meters)
        self.set_Moving_Window(1.5, 1, 512, 96) # in terms of plasma wavelength lambda_p
        # Simulation time
        self.set_sim_length(2e-3)
        # Plasma parameters
        self.set_plasma_size(1, 2, 2, 4)
        # Beam parameters
        self.set_beam_charge(250) # in pC
        # Moving window and diagnostic
        self.v_window = c
        # Initialize Plasma
        self.initialise_plasma()
        # Initialize Moving Window
        self.sim.set_moving_window(v=self.v_window)

    # --------------------
    # Simulation parameters
    # --------------------
    def set_Sim_control(self):
        self.use_cuda = True
        self.n_order = -1  # Infinite-order solver
        self.Nm = 2  # Azimuthal modes

    # --------------------
    # Calculating lambda_p
    # --------------------
    def get_lambda_p(self):
        self.omega_p = np.sqrt(self.n0 * e ** 2 / (m_e * epsilon_0))
        self.lambda_p = 2 * np.pi * c / self.omega_p
        print(f"lambda_p: {self.lambda_p*1e6} micron")

    # --------------------
    # Moving window domain size (meters)
    # --------------------
    def set_Moving_Window(self, zmax, rmax, Nz, Nr):
        self.zmax = zmax * self.lambda_p  # Moving window length
        self.zmin = -zmax * self.lambda_p  # Start window centered around z=0
        self.rmax = rmax * self.lambda_p  # Radial extent = 1 plasma wavelength
        self.Nz = Nz
        self.Nr = Nr
        self.dz = (self.zmax - self.zmin) / self.Nz
        self.dr = self.rmax / self.Nr
        # More conservative time step
        self.dt = 0.5 * self.dz / c  # CFL condition for electromagnetic codes
        print(f"Moving window size {2*zmax*1e6}µm *{2*rmax*1e6}µm with {Nz} * {Nr} grids")
    # --------------------
    # Simulation time
    # --------------------
    def set_sim_length(self, total_propagation):
        self.total_propagation = total_propagation  # 20 cm
        self.t_final = self.total_propagation / c  # Time for window to travel 20 cm
        self.N_steps = int(self.t_final / self.dt)
        print(f"Simulation steps: {self.N_steps}")

    # --------------------
    # Plasma parameters
    # --------------------
    def set_plasma_size(self, p_zmin, p_nz, p_nr, p_nt):
        self.p_zmin = p_zmin * self.lambda_p
        self.p_zmax = self.total_propagation  # Fill the moving window
        self.p_rmax = self.rmax  # Fill radial extent (1*lambda_p)
        self.p_nz = p_nz  # Further reduced particles
        self.p_nr = p_nr
        self.p_nt = p_nt
        print(f"Simulation start: {self.p_zmin}, with grids z: {self.p_nz}, r: {self.p_nr}, t: {self.p_nt}")

    # --------------------
    # Beam parameters
    # --------------------
    def set_beam_charge(self, beam_charge):
        self.total_charge = beam_charge  # in pC
        print(f"Beam charge: {self.total_charge} pC")

    # --------------------
    # Initialize simulation
    # --------------------
    def initialise_plasma(self):
        print("Initializing simulation(1)...")
        self.sim = Simulation(self.Nz, self.zmax, self.Nr, self.rmax, self.Nm, self.dt, zmin=self.zmin,
                              n_order=self.n_order, use_cuda=self.use_cuda,
                              boundaries={'z': 'open', 'r': 'reflective'}, )

        self.electrons = self.sim.add_new_species(
            q=-e, m=m_e, n=self.n0,
            dens_func=None,
            p_zmin=self.p_zmin,  # Start where plasma actually exists
            p_zmax=self.p_zmax,
            p_rmax=self.p_rmax,
            p_nz=self.p_nz, p_nr=self.p_nr, p_nt=self.p_nt
        )

        # Plasma ions
        self.protons = self.sim.add_new_species(
            q=e, m=m_p, n=self.n0,
            dens_func=None,
            p_zmin=self.p_zmin,  # Match electrons
            p_zmax=self.p_zmax,
            p_rmax=self.p_rmax,
            p_nz=self.p_nz, p_nr=self.p_nr, p_nt=self.p_nt
        )


    def set_input_Gaussian(self, sigma_z =1.4e-5, sigma_r =1.85e-5, n_emit=4.426719e-6, sig_gamma=6.05, n_macroparticles=262144, injection_plane= 1):
        self.gamma_b = self.total_charge / 0.511  # Mev/0.511 Reduced from extremely high value
        # Make beam much more compact: ~10 µm instead of ~84 µm
        self.beam_sigma_z = sigma_z  # Longitudinal RMS size ~8 µm
        self.beam_sigma_r = sigma_r  # Transverse RMS size ~8 µm
        self.beam_z0 = 1/2 * self.lambda_p    # Start at 105.6 µm (before plasma)
        self.n_particles = n_macroparticles

        self.beam = add_particle_bunch_gaussian(
                self.sim, q=-e, m=m_e, gamma0=self.gamma_b,
                n_emit=n_emit,  # Zero emittance (idealized beam)
                sig_r=self.beam_sigma_r, sig_z=self.beam_sigma_z, sig_gamma=sig_gamma,
                n_physical_particles= int(self.total_charge*1e-12/e),
                n_macroparticles= n_macroparticles,
                tf=0,  # Injection time
                zf=self.beam_z0,  # Beam center position
                z_injection_plane= injection_plane * self.lambda_p,
                boost=None,)
        print(f"Simulation Build")
        print(f"Number of MacroParticle: {self.n_particles} ")
        print(f"gamma: {self.gamma_b} ")
        print(f"Longitudinal RMS(sigma_z): {self.beam_sigma_z*1e6} µm")
        print(f"Transverse RMS(sigma_r): {self.beam_sigma_r*1e6} µm")
        print(f"normalised emittance: {n_emit*1e6}µm*mrad")

    
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
        print(f"weight = {self.w_beam}")
        self.beam = add_particle_bunch_from_arrays(self.sim, q=-e, m=m_e, x=self.x_beam,
                                                    y=self.y_beam, z=self.z_beam, ux=self.px_beam, 
                                                    uy =self.py_beam, uz=self.pz_beam, w=self.w_beam,
                                                    z_injection_plane = 1 * self.lambda_p)
        

    def run(self, outputdir, diag_period = 0, fieldtype =["E", "rho"], extra_species = None):
        self.write_dir = outputdir

        if diag_period > 0:
            self.diag_period = diag_period
        else:
            self.diag_period = self.N_steps

        species = {"electrons": self.electrons}
        if extra_species is not None:
            species.update(extra_species)
        # --------------------
        # Diagnostics
        # -------------------
        self.sim.diags = [
            # Save only E-field components (not B-field)
            FieldDiagnostic(self.diag_period, self.sim.fld, comm=self.sim.comm,
                fieldtypes=fieldtype,  # Only E-field and charge density
                write_dir=self.write_dir,

            ),
            
            # Save all particle species
            ParticleDiagnostic(self.diag_period, species, comm=self.sim.comm,
                write_dir=self.write_dir,
            )
        ]
        # --------------------
        # Run simulation
        # --------------------
        print("Starting simulation...")
        print(f"diag_period: {self.diag_period} ")
        print(f"Field type: {fieldtype}")
        print(f"species: {species}")
        try:
            self.sim.step(self.N_steps+1)
            print("Simulation completed successfully!")
        except Exception as e:
            print(f"Simulation failed with error: {e}")
            print("Try reducing the time step or beam density further.")
                
    def get_output_sdds(self,outputdir,species):
        input = f"{self.write_dir}/hdf5/data{self.diag_period:08d}.h5"
        pyclara.Converters.fbpic2sdds(input, outputdir, species)




