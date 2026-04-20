import numpy as np
import pyclara
import h5py
from scipy.constants import c, e, m_e, m_p, epsilon_0
from fbpic.main import Simulation
from fbpic.openpmd_diag import FieldDiagnostic, ParticleDiagnostic
from fbpic.lpa_utils.bunch import add_particle_bunch_gaussian
from fbpic.lpa_utils.bunch import add_particle_bunch_from_arrays



class Fbpic_runner:
    """
    Runs FBPIC plasma wakefield simulations as part of an S2E pipeline.

    Sets up and runs a particle-in-cell (PIC) simulation of a plasma
    wakefield accelerator using FBPIC. Handles plasma initialisation,
    beam injection, running the simulation, and returning output particles
    in OpenPMD HDF5 format.

    Parameters
    ----------
    outputdir : str
        Path to the directory where simulation output HDF5 files will
        be written.

    Attributes
    ----------
    n0 : float
        Background plasma density [m^-3].
    lambda_p : float
        Plasma wavelength [m].
    omega_p : float
        Plasma frequency [rad/s].
    total_propagation : float
        Total propagation distance of the simulation [m].
    N_steps : int
        Total number of simulation time steps.
    write_dir : str
        Directory where output HDF5 files are written.
    sim : fbpic.Simulation
        The FBPIC simulation object.
    electrons : fbpic species
        Background plasma electron species.
    protons : fbpic species
        Background plasma ion (proton) species.
    beam : fbpic species
        Beam (driver) particle species.
    Examples
    --------
    # f = Fbpic_runner('/path/to/output')
    # f.set_plasma_density(1e22)
    # output = f.track(elegant_output_h5)
    """
    def __init__(self, outputdir):
        # Simulation parameters
        self.set_Sim_control()
        # Plasma parameters
        self.set_plasma_density(1e22)  # Background plasma density [m^-3] (example value)
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
        self.write_dir = outputdir
        self.initialise_plasma()
        # Initialize Moving Window
        self.sim.set_moving_window(v=self.v_window)

    # --------------------
    # Simulation parameters
    # --------------------
    def set_Sim_control(self, n_order = -1, Nm = 2 ):
        """
                Set simulation control parameters.

                Sets CUDA usage, spectral solver order, and number of
                azimuthal modes.
        """
        self.use_cuda = True
        self.n_order = n_order  # Infinite-order solver
        self.Nm = Nm  # Azimuthal modes

    # --------------------
    # Set plasma density
    # --------------------
    def set_plasma_density(self, n0):
        """
       Set the background plasma density.

       Parameters
       ----------
       n0 : float
           Background plasma electron density [m^-3].
        """
        self.n0 = n0  # Background plasma density [m^-3] (example value)
        print(f"Plasma density: {self.n0} [m^-3]")

    # --------------------
    # Calculating lambda_p
    # --------------------
    def get_lambda_p(self):
        """
        Calculate and store the plasma wavelength and frequency.

        Computes omega_p and lambda_p from the plasma density n0.
        Must be called after set_plasma_density().
        """
        self.omega_p = np.sqrt(self.n0 * e ** 2 / (m_e * epsilon_0))
        self.lambda_p = 2 * np.pi * c / self.omega_p
        print(f"lambda_p: {self.lambda_p*1e6} micron")

    # --------------------
    # Moving window domain size (meters)
    # --------------------
    def set_Moving_Window(self, zmax, rmax, Nz, Nr):
        """
        Set the moving window size and grid resolution.

        All spatial arguments are given in units of the plasma
        wavelength lambda_p.

        Parameters
        ----------
        zmax : float
            Half-length of the moving window in units of lambda_p.
            Window spans [-zmax, +zmax] * lambda_p.
        rmax : float
            Radial extent of the window in units of lambda_p.
        Nz : int
            Number of grid cells in the longitudinal direction.
        Nr : int
            Number of grid cells in the radial direction.
        """
        self.zmax = zmax * self.lambda_p  # Moving window length
        self.zmin = -zmax * self.lambda_p  # Start window centered around z=0
        self.rmax = rmax * self.lambda_p  # Radial extent = 1 plasma wavelength
        self.Nz = Nz
        self.Nr = Nr
        self.dz = (self.zmax - self.zmin) / self.Nz
        self.dr = self.rmax / self.Nr
        # More conservative time step
        self.dt = 0.5 * self.dz / c  # CFL condition for electromagnetic codes
        print(f"Moving window size {2*self.zmax*1e6}µm *{2*self.rmax*1e6}µm with {Nz} * {Nr} grids")
    # --------------------
    # Simulation time
    # --------------------
    def set_sim_length(self, total_propagation):
        """
                Set the total propagation distance and calculate simulation steps.

                Parameters
                ----------
                total_propagation : float
                    Total distance the beam propagates through plasma [m].

                N_steps : int
                    Total propagation steps.
        """
        self.total_propagation = total_propagation  # 20 cm
        self.t_final = self.total_propagation / c  # Time for window to travel 20 cm
        self.N_steps = int(self.t_final / self.dt)
        print(f"Simulation steps: {self.N_steps}")

    # --------------------
    # Plasma parameters
    # --------------------
    def set_plasma_size(self, p_zmin, p_nz, p_nr, p_nt):
        """
        Set the plasma extent and macro-particle grid resolution.

        Parameters
        ----------
        p_zmin : float
            Longitudinal start of the plasma in units of lambda_p.
        p_nz : int
            Number of macro-particles per cell in the z direction.
        p_nr : int
            Number of macro-particles per cell in the r direction.
        p_nt : int
            Number of macro-particles per cell in the theta direction.
        """
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
        """
        Set the total beam charge.

        Parameters
        ----------
        beam_charge : float
            Total beam charge [pC].
        """
        self.total_charge = beam_charge  # in pC
        print(f"Beam charge: {self.total_charge} pC")

    # --------------------
    # Initialize simulation
    # --------------------
    def initialise_plasma(self):
        """
        Initialise the FBPIC simulation and add plasma species.

        Creates the FBPIC Simulation object and adds background
        electron and proton (ion) species filling the plasma volume.
        Must be called after all simulation parameters are set.
        """
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
        """
        Inject a Gaussian witness beam into the simulation.

        Parameters
        ----------
        sigma_z : float, optional
            Longitudinal RMS bunch length [m]. Default is 1.4e-5.
        sigma_r : float, optional
            Transverse RMS beam size [m]. Default is 1.85e-5.
        n_emit : float, optional
            Normalised transverse emittance [m·rad]. Default is 4.426719e-6.
        sig_gamma : float, optional
            RMS energy spread in units of gamma. Default is 6.05.
        n_macroparticles : int, optional
            Number of macro-particles. Default is 262144.
        injection_plane : float, optional
            z position of the injection plane in units of lambda_p.
            Default is 1.
        """
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
        """
        Inject a beam from a dictionary of particle arrays.

        Reads particle positions and momenta from a dictionary,
        applies a z-offset to centre the beam at half a plasma
        wavelength, and adds the beam species to the simulation.

        Parameters
        ----------
        input : dict
            Dictionary with keys 'x', 'y', 'z', 'px', 'py', 'p',
            each containing a numpy array of particle coordinates
            or momenta. Positions in [m], momenta as βγ.
        """
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

    def set_input_sdds(self, sdds_dir):
        """
        Inject a beam from an Elegant SDDS output file.

        Converts the SDDS file to a particle dictionary using
        sdds2fbpic and passes it to set_input_dict().

        Parameters
        ----------
        sdds_dir : str
            Path to the Elegant output SDDS file.
        """
        f = pyclara.Converters.sdds2fbpic(sdds_dir)
        self.set_input_dict(f)

    def set_input_h5(self, input):
        """
        Inject a beam from an OpenPMD HDF5 file.

        Intended to receive the output of Elegant_runner.get_output()
        and pass it to set_input_dict().

        Parameters
        ----------
        input : dict
            Particle dictionary in the format returned by sdds2fbpic,
            with keys 'x', 'y', 'z', 'px', 'py', 'p'.
        """
        self.set_input_dict(input)

    def run(self, diag_period = 0, fieldtype =["E", "rho"], extra_species = None):
        """
        Run the FBPIC simulation with diagnostics.

        Sets up field and particle diagnostics, then runs the simulation
        for N_steps time steps.

        Parameters
        ----------
        diag_period : int, optional
            Number of steps between diagnostic outputs. If 0 (default),
            only the final step is saved.
        fieldtype : list of str, optional
            Field components to save. Default is ['E', 'rho'].
        extra_species : dict, optional
            Additional particle species to include in diagnostics,
            as {name: species} pairs. Default is None.
        """
        if diag_period > 0:
            self.diag_period = diag_period
        else:
            self.diag_period = self.N_steps

        species = {"electrons": self.electrons, "beam": self.beam}
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
                
    def get_output_h5(self):
        """
        Read the final output HDF5 file and correct the z-offset.

        Opens the output HDF5 file from the last diagnostic step,
        restores the original z coordinates by reversing the offset
        applied in set_input_dict(), and returns the open file object.

        Returns
        -------
        h5py.File
            Output particle distribution in OpenPMD HDF5 format with
            z positions corrected back to the lab frame.
        """
        output = f"{self.write_dir}/hdf5/data{self.diag_period:08d}.h5"
        data = h5py.File(output, "r+")
        iteration = list(data['data'].keys())[0]  # get last iteration
        for species in data[f'data/{iteration}/particles'].keys():
            path = f'data/{iteration}/particles/{species}/position/z'
            data[path][:] = data[path][:] + self.beam_z0_off
        return data


    def track(self, particles):
        """
        Run a full FBPIC tracking step in the S2E pipeline.

        Convenience method that calls set_input_h5(), run(), and
        get_output_h5() in sequence. This is the method called by
        the S2E pipeline.

        Parameters
        ----------
        particles : dict
            Input particle dictionary in the format returned by
            Elegant_runner.get_output(), with keys 'x', 'y', 'z',
            'px', 'py', 'p'.

        Returns
        -------
        h5py.File
            Output particle distribution in OpenPMD HDF5 format, ready
            to be passed to the next tracker in the pipeline.

        Examples
        --------
        # output = f.track(elegant_output)
        """
        self.set_input_h5(particles)
        self.run()
        output = self.get_output_h5()
        return output


