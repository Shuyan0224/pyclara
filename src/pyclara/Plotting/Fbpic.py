import matplotlib.pyplot as plt
import numpy as np 
import h5py
from scipy.constants import c, e, m_e

def All(hdf5file,species): 

    plt.figure()
    
    plt.subplot(2, 3, 1)
    X(hdf5file,species)

    plt.subplot(2, 3, 2)
    Y(hdf5file,species)

    plt.subplot(2, 3, 3)
    Z(hdf5file,species)

    plt.subplot(2, 3, 4)
    Px(hdf5file,species)

    plt.subplot(2, 3, 5)
    Py(hdf5file,species)

    plt.subplot(2, 3, 6)
    Pz(hdf5file,species)

    plt.tight_layout()

def X(hdf5file,species):
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        x = f[f"/data/{iteration}/particles/{species}/position/x"][:]
        plt.hist(x * 1e3, bins=100)
        plt.xlabel("x position (mm)")
        plt.ylabel("No. of macroparticles")

def Y(hdf5file,species):
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        y = f[f"/data/{iteration}/particles/{species}/position/y"][:]
        plt.hist(y * 1e3, bins=100)
        plt.xlabel("y position (mm)")
        plt.ylabel("No. of macroparticles")

def Z(hdf5file,species):
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        z = f[f"/data/{iteration}/particles/{species}/position/z"][:]
        plt.hist((z - np.mean(z)) * 1e3, bins=100)
        plt.xlabel("z position (mm)")
        plt.ylabel("No. of macroparticles")

def Px(hdf5file,species):
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        px = f[f"/data/{iteration}/particles/{species}/momentum/x"][:]
        plt.hist(px, bins=100)
        plt.yscale("log")
        plt.xlabel("$u_x$ (β·γ)")
        plt.ylabel("No. of macroparticles")

def Py(hdf5file,species):
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        py = f[f"/data/{iteration}/particles/{species}/momentum/y"][:]
        plt.hist(py, bins=100)
        plt.yscale("log")
        plt.xlabel("$u_y$ (β·γ)")
        plt.ylabel("No. of macroparticles")

def Pz(hdf5file,species):
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        pz = f[f"/data/{iteration}/particles/{species}/momentum/z"][:]
        plt.hist(pz, bins=100)
        plt.yscale("log")
        plt.xlabel("$u_z$ (β·γ)")
        plt.ylabel("No. of macroparticles")

def Wakefield_Ez(hdf5file,species):
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        attrs = f[f'data/{iteration}/fields/E'].attrs
        dr, dz = attrs['gridSpacing']
        r0, z0 = attrs['gridGlobalOffset']
        Ez = f[f"/data/{iteration}/fields/E/z"][:]
        Ez = Ez[0] + 2 * Ez[1]
        Nr, Nz = Ez.shape
        z = (z0 + np.arange(Nz) * dz) * 1e6  # µm
        r = (r0 + np.arange(Nr) * dr) * 1e6  # µm
        Ez_full = np.vstack([Ez[::-1], Ez])
        r_full = np.concatenate([-r[::-1], r])
        v = np.percentile(np.abs(Ez_full / 1e9), 99)
        plt.figure(figsize=(10, 4))
        scale = np.max(np.abs(Ez[0])) / (r_full[-1] * 0.8)
        plt.plot(z, Ez[0, :] / scale, color='black', lw=1.5)
        plt.imshow(Ez_full / 1e9, origin='lower', aspect='auto', cmap='RdBu_r',
                   vmin=-v, vmax=v, extent=[z[0], z[-1], r_full[0], r_full[-1]])
        plt.colorbar(label='GV/m', pad=0.10)
        plt.xlabel('z (µm)');
        plt.ylabel('r (µm)');
        plt.title('longitudinal wakefield')
        original_ylim = plt.ylim()
        ax2 = plt.twinx()
        ax2.set_ylim(np.array(original_ylim) * scale / 1e9)
        ax2.set_ylabel('Ez at r=0 (GV/m)')
        ax2.tick_params(axis='y')
        plt.tight_layout()
        plt.show()

def rho (hdf5file,species):
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        attrs = f[f'data/{iteration}/fields/E'].attrs
        dr, dz = attrs['gridSpacing']
        r0, z0 = attrs['gridGlobalOffset']
        w = f[f"/data/{iteration}/particles/electrons/weighting"][:]
        mean_w = np.mean(w)
        rho = f[f'data/{iteration}/fields/rho'][:]
        rho = rho[0] + 2 * rho[1]
        Nr, Nz = rho.shape
        z = (z0 + np.arange(Nz) * dz) * 1e6
        r = (r0 + np.arange(Nr) * dr) * 1e6
        rho_full = np.vstack([rho[::-1], rho]) / -e / 1e6 * mean_w
        r_full = np.concatenate([-r[::-1], r])
        v = np.percentile(np.abs(rho_full), 99)
        plt.figure(figsize=(10, 4))
        plt.imshow(rho_full, origin='lower', aspect='auto', cmap='Blues',
                   vmin=-v, vmax=v, extent=[z[0], z[-1], r_full[0], r_full[-1]])
        plt.colorbar(label='C/cm³')
        plt.xlabel('z (µm)');
        plt.ylabel('r (µm)');
        plt.title('The electron density')
        plt.tight_layout()
        plt.show()

def phase_space(hdf5file,species):
    MeV = 1e6 * e
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        beam = f[f'data/{iteration}/particles/{species}']
        px = beam['momentum/x'][:]
        py = beam['momentum/y'][:]
        pz = beam['momentum/z'][:]
        x = beam['position/x'][:] * 1e6
        y = beam['position/y'][:] * 1e6
        z = beam['position/z'][:] * 1e6
        w = beam['weighting'][:]

    p_norm = np.sqrt(px ** 2 + py ** 2 + pz ** 2) / (m_e * c)
    gamma = np.sqrt(1 + p_norm ** 2)
    energy = (gamma - 1) * m_e * c ** 2 / MeV

    px_norm = px / (m_e * c)
    py_norm = py / (m_e * c)

    # --- x-px ---
    plt.figure(figsize=(6, 4))
    plt.scatter(x, px_norm, s=0.5, c='steelblue', alpha=0.3)
    plt.xlabel('x (µm)');
    plt.ylabel('$p_x / m_e c$')
    plt.title('x-px phase space')
    plt.tight_layout()

    # --- y-py ---
    plt.figure(figsize=(6, 4))
    plt.scatter(y, py_norm, s=0.5, c='steelblue', alpha=0.3)
    plt.xlabel('y (µm)');
    plt.ylabel('$p_y / m_e c$')
    plt.title('y-py phase space')
    plt.tight_layout()

    # --- z-energy ---
    plt.figure(figsize=(6, 4))
    plt.scatter(z, energy, s=0.5, c='steelblue', alpha=0.3)
    plt.xlabel('z (µm)');
    plt.ylabel('Energy (MeV)')
    plt.title('z-Energy phase space')
    plt.tight_layout()

    plt.show()

def divergence(hdf5file,species):
    with h5py.File(hdf5file, 'r') as f:
        iteration = list(f["data"].keys())[0]
        beam = f[f'data/{iteration}/particles/{species}']
        px = beam['momentum/x'][:]
        py = beam['momentum/y'][:]
        pz = beam['momentum/z'][:]
        x = beam['position/x'][:] * 1e6
        y = beam['position/y'][:] * 1e6
        xp = px / pz * 1e3  # mrad
        yp = py / pz * 1e3  # mrad
        # --- x-x' ---
        plt.figure(figsize=(6, 4))
        plt.scatter(x, xp, s=0.5, c='steelblue', alpha=0.3)
        plt.xlabel('x (µm)');
        plt.ylabel("x' (mrad)")
        plt.title("x-x' divergence")
        plt.tight_layout()

        # --- y-y' ---
        plt.figure(figsize=(6, 4))
        plt.scatter(y, yp, s=0.5, c='steelblue', alpha=0.3)
        plt.xlabel('y (µm)');
        plt.ylabel("y' (mrad)")
        plt.title("y-y' divergence")
        plt.tight_layout()

        plt.show()