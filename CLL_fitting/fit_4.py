import numpy as np
import matplotlib.pyplot as plt

def main():
    global sigma, Mmass, E_in

    # ==============================================================
    # RNG + pre-generated samples
    # ==============================================================
    rng = np.random.default_rng(42)
    X_fixed = rng.random((200000, 6))

    # ==============================================================
    # Parameters
    # ==============================================================
    Tw = 300
    M = 1
    kb = 1.380649E-23
    m = 2.65E-26 #AO
    # m = 6.63E-26 #Argon
    Mmass = 0.015999
    Rg = kb / m
    # Rg = 8314.5

    V = 7930 #8000#7888
    # V = 900
    theta_beam = 70  # degrees

    # E_in = 0.5 * Mmass * V**2
    E_in = 0.5 * Mmass * V**2

    FWHM = 614.44744
    sigma = 0.7272475381313256 * FWHM

    # ==============================================================
    # Load data
    # ==============================================================
    exp_E_in = 503.1 #495.84234
    data = np.loadtxt(f"Austin_Scattering_1_{theta_beam}_incedent_angle.csv",
                      delimiter=",", skiprows=1)
    exp_angles = data[:, 0]
    exp_count = data[:, 3]
    exp_energy = data[:, 1] / exp_E_in
    exp_count /= np.max(exp_count)

    # exp_E_in = 495.84234
    # data = np.loadtxt(f"coated.csv",
    #                   delimiter=",", skiprows=1)
    # exp_angles = data[:, 0]
    # exp_count = data[:, 2]
    # exp_energy = data[:, 1] / exp_E_in
    # exp_count /= np.max(exp_count)

    exp2 = 0

    if exp2 == 1:
        exp_E_in = 495.84234
        data = np.loadtxt(f"uncoated.csv",
                        delimiter=",", skiprows=1)
        exp_angles_2 = data[:, 0]
        exp_count_2 = data[:, 2]
        exp_energy = data[:, 1] / exp_E_in
        exp_count_2 /= np.max(exp_count_2)

    new_angles = np.linspace(-90, 90, 181)

    # ==============================================================
    # Model
    # ==============================================================
    def compute_distribution(alpha_n, alpha_t, Tw, M):

        # -------- velocities --------
        V_thermal_n = rng.normal(0.0, sigma, X_fixed.shape[0])
        V_thermal_t1 = rng.normal(0.0, sigma, X_fixed.shape[0])

        xi_t1_i = V * np.sin(np.radians(theta_beam)) + V_thermal_t1
        xi_n_i  = -V * np.cos(np.radians(theta_beam)) + V_thermal_n

        xi_mp_w = np.sqrt((2 * Rg * Tw) / M)
        # xi_mp_w = np.sqrt((1039*Tw))

        # print(2 * Rg / M)
        # print(Rg, Tw, M)

        r1 = np.sqrt(-alpha_n * np.log(X_fixed[:, 0]))
        phi2 = 2 * np.pi * X_fixed[:, 1]
        xi_n_m = (xi_n_i / xi_mp_w) * np.sqrt(1 - alpha_n)
        xi_n_r = xi_mp_w * np.sqrt(r1**2 + xi_n_m**2 + 2*r1*xi_n_m*np.cos(phi2))

        r3 = np.sqrt(-alpha_t * np.log(X_fixed[:, 2]))
        phi4 = 2 * np.pi * X_fixed[:, 3]
        xi_t1_m = (xi_t1_i / xi_mp_w) * np.sqrt(1 - alpha_t)
        xi_t1_r = xi_mp_w * (xi_t1_m + r3 * np.cos(phi4))

        r5 = np.sqrt(-alpha_t * np.log(X_fixed[:, 4]))
        phi6 = 2 * np.pi * X_fixed[:, 5]
        xi_t2_r = xi_mp_w * r5 * np.cos(phi6)

        vx, vy, vz = xi_t1_r, xi_n_r, xi_t2_r
        vel = np.column_stack((vx, vy, vz))

        # -------- angular distribution --------
        norm = np.array([0.0, 1.0, 0.0])
        tang = np.array([1.0, 0.0, 0.0])

        theta_mid, counts, bins, scatter_angle = angular_scattering_profile(
            vel,
            norm,
            tangent=tang,
            nbins=180,
            degrees=True
        )

        # Interpolate count
        count = np.interp(new_angles, theta_mid, counts, left=0, right=0)
        count /= np.max(count)

        # -------- energy ratio --------
        # scatter_angle is radians → convert to degrees
        scatter_angle_deg = np.degrees(scatter_angle)

        # bins must match angular_scattering_profile
        bin_edges_deg = np.linspace(-90, 90, len(theta_mid) + 1)

        bin_indices = np.digitize(scatter_angle_deg, bin_edges_deg)

        energy_ratio = np.zeros(len(theta_mid))

        for i in range(len(theta_mid)):
            idx = np.where(bin_indices == i + 1)[0]
            if idx.size == 0:
                continue

            # v2 = vx[idx]**2 + vy[idx]**2 + vz[idx]**2
            v2 = vx[idx]**2 + vy[idx]**2
            E_out = 0.5 * Mmass * v2
            energy_ratio[i] = np.mean(E_out / E_in)

        energy_ratio = np.interp(new_angles, theta_mid, energy_ratio, left=0, right=0)

        return count, energy_ratio

    # ==============================================================
    # Runs
    # ==============================================================
    # alpha_n = np.array([0.3, 0.2, 0.2, 0.2, 1])
    # alpha_t = np.array([0.7, 0.8, 0.8, 0.8, 1])
    # wall_temp = np.array([300*50, 300*50, 300*60, 300*40, 300*50])
    # Mvec = np.array([1, 1, 1, 1, 1])

    # alpha_n = np.array([0.3, 0.05, 0.3, 1, 0, 1])
    # alpha_t = np.array([0.7, 0.9, 0.7, 0.7, 0, 1])
    alpha_n = np.array([0.05, 0.1, 0.2, 0.3, 0, 1])
    alpha_t = np.array([0.9, 0.9, 0.9, 0.9, 0, 1])
    wall_temp = np.array([300*120, 300*120, 300*100, 300*100, 1, 1, 1])
    Mvec = np.array([1, 1, 1, 1, 1, 1])

    # alpha_n = np.array([0.3, 0.2, 0.25, 0, 1])
    # alpha_t = np.array([0.1, 0.7, 0.5, 0, 1])
    # wall_temp = np.array([1354.15, 1354.15, 1354.15, 1, 1, 1])
    # Mvec = np.array([1, 1, 1, 1, 1])


    # alpha_n = np.array([0.3, 0.18, 0.2, 0, 1])
    # alpha_t = np.array([0.58, 0.4, 0.95, 0, 1])
    # # Either
    # wall_temp = np.array([300/0.17, 300/0.15, 300/0.14, 1, 1, 1])
    # Mvec = np.array([1, 1, 1, 1, 1])
    # alpha_n = np.array([0.2])
    # alpha_t = np.array([0.95])
    # # Either
    # wall_temp = np.array([300/0.14])
    # Mvec = np.array([1])


    # Or 
    # wall_temp = np.array([300, 300, 300, 300, 300, 300, 300])   
    # Mvec = np.array([0.17, 0.15, 0.14, 1, 1])

    # wall_temp = np.array([300*50, 300*50, 300*50, 1, 1, 1])
    # Mvec = np.array([1, 1, 1, 1, 1])

    

    n = len(new_angles)
    m = len(alpha_n)

    fit = np.zeros((n, m))
    efit = np.zeros((n, m))

    for i in range(m):
        fit[:, i], efit[:, i] = compute_distribution(
            alpha_n[i], alpha_t[i], wall_temp[i], Mvec[i]
        )

    # ==============================================================
    # Plots
    # ==============================================================

    # Cartesian
    plt.figure(figsize=(7,5))
    # plt.plot(exp_angles, exp_count, "o", color='black', label="Experiment")
    plt.plot(exp_angles, exp_count, "o", color='black', label="Coated")
    if exp2 == 1:
        plt.plot(exp_angles_2, exp_count_2, "s", color='black', label="Uncoated")

    for i in range(m):
        plt.plot(new_angles, fit[:, i],label=f"a_n={alpha_n[i]:.2f}, a_t={alpha_t[i]:.2f}")

    plt.xlabel("Scattering angle (deg)")
    plt.ylabel("Normalized intensity")
    plt.legend()
    plt.savefig("model_cartesian.png", dpi=300)

    # Energy
    plt.figure(figsize=(5,5))
    plt.plot(exp_angles, exp_energy, "o", color='black')

    for i in range(m):
        plt.plot(new_angles, efit[:, i],label=f"a_n={alpha_n[i]:.2f}, a_t={alpha_t[i]:.2f}")
        

    plt.xlabel("Scattering angle (deg)")
    plt.ylabel("Energy ratio")
    plt.xlim(0, 90)
    plt.legend()
    plt.savefig("model_energy.png", dpi=300)

    # Polar
    plt.figure(figsize=(7,5))
    ax = plt.subplot(projection="polar")

    # ax.plot(np.radians(exp_angles), exp_count, "o", color='black')
    ax.plot(np.radians(exp_angles), exp_count, "o", color='black',markerfacecolor='none', label="Coated")
    if exp2 == 1:
        ax.plot(np.radians(exp_angles_2), exp_count_2, "s", color='black',markerfacecolor='none', label="Uncoated")


    ax.plot(np.radians(new_angles),
            np.cos(np.radians(new_angles)), '--', color='k')

    for i in range(m):
        ax.plot(np.radians(new_angles), fit[:, i],label=f"a_n={alpha_n[i]:.2f}, a_t={alpha_t[i]:.2f}")
   

    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    angles = np.arange(-90, 91, 10)
    ax.set_thetagrids(angles)
    plt.legend()
    plt.savefig("model_polar.png", dpi=300)

    # ==============================================================
    # Grid
    # ==============================================================
    # grid_n = 25
    # grid_t = 20
    # Tw = 300*50
    # alphas_n = np.linspace(0,1,grid_n)
    # alphas_t = np.linspace(0,1,grid_t)
    # err_surface = np.zeros((grid_n, grid_t))

    # for i, an in enumerate(alphas_n):
    #     for j, at in enumerate(alphas_t):
    #         gridfit, gridefit = compute_distribution(
    #             an, at, Tw, 1
    #             )

    #         gridfit = gridfit / np.max(gridfit)
    #         gridfit_int = np.interp(exp_angles, new_angles, gridfit)
    #         err_surface[i,j] = np.sqrt( np.mean((exp_count - gridfit_int)**2) )

    # # Plot RMS error surface
    # AN, AT = np.meshgrid(alphas_n, alphas_t, indexing='ij')

    # plt.figure(figsize=(6,5))
    # cs = plt.contourf(AN, AT, err_surface, levels=30, cmap='viridis')
    # # cs = plt.pcolormesh(AN, AT, err_surface[:-1, :-1], shading="flat", cmap="viridis")
    # cbar = plt.colorbar(cs)
    # cbar.set_label("RMS error")
    # plt.xlabel("α_n (normal)")
    # plt.ylabel("α_t (tangential)")
    # # Best fit marker
    # # plt.plot(alpha_n_fit, alpha_t_fit, 'ro', label="Best fit (0.3, 0.58)")
    # # Best found in grid
    # imin, jmin = np.unravel_index(np.argmin(err_surface), err_surface.shape)
    # grid_alpha_n = alphas_n[imin]
    # grid_alpha_t = alphas_t[jmin]
    # plt.plot(grid_alpha_n, grid_alpha_t, 'wx', ms=10, mew=2, label="Grid min")

    # plt.savefig("grid.png", dpi=300)


# ==============================================================
# DO NOT MODIFY
# ==============================================================

def angular_scattering_profile(
    V,
    normal,
    tangent=None,
    nbins=180,
    degrees=True,
    eps=1e-8
):
    V = np.asarray(V)
    n = np.asarray(normal)
    n = n / np.linalg.norm(n)

    if tangent is None:
        ref = np.array([1.0, 0.0, 0.0])
        if np.abs(np.dot(ref, n)) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        tangent = np.cross(n, ref)

    t = tangent / np.linalg.norm(tangent)

    Vn = V / np.linalg.norm(V, axis=1, keepdims=True)

    # cos_theta = np.clip(np.dot(Vn, n), -1.0, 1.0)
    # theta = np.arccos(cos_theta)

    theta = np.arctan2(V[:,0],V[:,1]) # what Minton did!
    # theta = np.arctan2(abs(V[:,0]),V[:,1])

    # print(np.min(theta),np.max(theta))

    phi = np.arctan2(V[:,2],V[:,0])
    mask1 = (phi < np.radians(5)) & (phi > np.radians(-5)) | (phi > np.radians(180-5)) & (phi < np.radians(180+5))
    # print(np.min(phi),np.max(phi))
    # theta = theta[mask1]

    sign = np.sign(np.dot(Vn, t))
    sign[sign == 0] = 1.0

    theta_signed = theta

    # print(np.nonzero(theta),np.nonzero(theta_signed))

    # theta_signed = theta * sign
    # mask = np.abs(theta_signed) <= np.pi / 2
    # theta_signed = theta_signed[mask]

    # print(np.min(theta_signed),np.max(theta_signed))

    bins = np.linspace(-np.pi/2, np.pi/2, nbins + 1)
    hist, edges = np.histogram(theta_signed, bins=bins)

    theta_centers = 0.5 * (edges[:-1] + edges[1:])
    profile = hist

    if degrees:
        theta_centers = np.degrees(theta_centers)

    return theta_centers, profile, bins, theta_signed


if __name__ == "__main__":
    main()