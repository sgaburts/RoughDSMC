
# Computes surface drag
# Computes scattering distribution from particles
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
import math

def main():
    cases = ["new_cll70_50","new_cll70_60","new_cll70_70","."]
    # cases = ["cll_50","cll_60_incomplete","cll_70","cll_80"]
    # cases = ["spec_70_0.1nm","spec_70_0.5nm","spec_70_1nm","spec_70_2nm","spec_70_3nm","spec_70_5nm","spec_70_7nm","spec_70_9nm","spec_70_10nm"]

    theta_beam = [50, 60, 70, 80, 70, 70, 70, 70, 70, 70]
    lines = ['-','-','-','-','-','-','-','-','-','-','-','-','-','-.']
    col = ['b','b','b','r','r','r','b','b','b','r','r','r','r','b','b','b']

    n = len(cases)

    count = [None]*n
    theta_centers = [None]*n

    fig1, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    for i in range(n):
        casename = cases[i]
        zlim = 10

        # ------------------------------------------------------------
        # Find all scatter files automatically
        pattern = os.path.join(casename, "scatter.*.dat")
        files = glob.glob(pattern)



        if len(files) == 0:
            print(f"No files found for {casename}")
            continue

        # Sort files by step number
        def extract_step(f):
            match = re.search(r"scatter\.(\d+)\.dat", f)
            return int(match.group(1)) if match else -1

        files = [f for f in files if extract_step(f) > 0]
        files = sorted(files, key=extract_step)

        # ------------------------------------------------------------
        count_sum = None
        theta_centers_sum = None
        C_sum = None
        E_sum = None
        

        valid_file_count = 0

        for f in files:
            try:
                c, tc, C2, E2, theta_mid, phi_bins = scattering_profile(zlim, f, n_theta=180, n_phi=360)
            except Exception as e:
                print(f"Skipping {f}: {e}")
                continue

            if count_sum is None:
                count_sum = np.array(c, dtype=float)
                theta_centers_sum = np.array(tc, dtype=float)
                C_sum = np.array(C2, dtype=float)
                E_sum = np.array(E2, dtype=float)
            else:
                count_sum += c
                theta_centers_sum += tc
                C_sum += C2
                E_sum += E2

            valid_file_count += 1

        # ------------------------------------------------------------
        if valid_file_count > 0:
            count[i] = count_sum / valid_file_count
            theta_centers[i] = theta_centers_sum / valid_file_count
            C_avg = C_sum / valid_file_count
            E_avg = E_sum / valid_file_count
        else:
            print(f"No valid files processed for {casename}")


        Fx, Fy, Fz, A, An, th = compute_drag(casename+"/prop_surfos.1000.dat")
        # Fx = 0
        # Fy = 0
        # Fz = 0
        # A = 0
        # An = 0
        # th = 0


        print(casename,'(Fx, Fy, Fz, A, A_n, th) = (',Fx,Fy,Fz,A,An,th,')')

        if i == 0:
            with open('scattering.txt', "a") as writefile:  # Open file in append mode
                line = "theta_beam\t"+"Fx\t"+"Fy\t"+"Fz\t"+ "\t".join(f"{x:.6f}" for x in np.degrees(theta_centers[i])) + "\n"
                writefile.write(line)
        with open('scattering.txt', "a") as writefile:  # Open file in append mode
            line = f"{theta_beam[i]}\t" + f"{Fx}\t" + f"{Fy}\t" + f"{Fz}\t" + "\t".join(f"{x:.6f}" for x in count[i]) + "\n"
            writefile.write(line)

        #------------------------------------------------------------
        # Experimental data
        # data = np.loadtxt("Austin_Scattering_1_70_incedent_angle.csv", delimiter=",", skiprows=1)
        # exp_angles = data[:,0]
        # exp_count = data[:,3]
        # # exp_energy = data[:,1]/exp_E_in
        # # Normalize exp_count
        # exp_count = exp_count/np.max(exp_count)


        #------------------------------------------------------------
        # Polar plot
        count[i] = count[i] / np.max(count[i])

        # if i==0:
        #     ax.plot(np.radians(exp_angles), exp_count, "o", color='black', label="Experiment")
        #     ax.plot(theta_centers[i], np.cos(theta_centers[i]), linestyle='--',color='k', label='cosine')

        ax.plot(theta_centers[i], count[i], linestyle=lines[i], label=casename)

        ax.set_thetamin(-90)
        ax.set_thetamax(90)
        ax.set_theta_zero_location("N")    # Optional: 0° at top
        ax.set_theta_direction(-1)         # Optional: clockwise
        angles = np.arange(-90, 90, 10)  # every 10 degrees
        ax.set_thetagrids(angles)
        plt.tight_layout()  # Adjust layout to fit everything
        # ax.legend(loc='center left', bbox_to_anchor=(0.8, 0.8))
        ax.legend()
        fig1.savefig("plot.png", dpi=300)

        # 2D plots
        # 2D plots

        C_avg /= np.max(C_avg)
        E_avg /= np.max(E_avg)

        # --- edges ---
        phi_edges = phi_bins                              # (361,)
        theta_edges = np.linspace(0, np.pi/2, C_avg.shape[1] + 1)  # (181,)
        theta_edges_deg = np.degrees(theta_edges)

        # --- wrap phi (periodic) ---
        C_plot = np.vstack([C_avg, C_avg[0:1, :]])        # (361, 180)

        # --- meshgrid (FIX ORDER HERE) ---
        theta_grid, phi_grid = np.meshgrid(theta_edges_deg, phi_edges)

        plt.figure(figsize=(7,5))
        ax2 = plt.subplot(projection="polar")

        cs = ax2.pcolormesh(phi_grid, theta_grid, C_plot[:-1,:], shading='auto')

        plt.colorbar(cs, label="Intensity")
        ax2.set_title(f"2D Scattering Intensity: {casename}")

        ax2.set_ylim(0, 90)

        plt.tight_layout()
        plt.savefig(f"{casename}_2D_intensity.png", dpi=300)

    


def scattering_profile(zlim, datapath, n_theta=180, n_phi=360):
    import numpy as np

    data = np.loadtxt(datapath, skiprows=9)
    vx, vy, vz = data[:, 5], data[:, 6], data[:, 7]
    x, y, z = data[:, 2], data[:, 3], data[:, 4]

    xplate = -2.4951e-06
    mask = x > xplate

    vx = vx[mask]
    vy = vy[mask]
    vz = vz[mask]

    vel = np.column_stack((vx, vy, vz))

    # ------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------
    norm = np.array([0.0, 1.0, 0.0])
    tang = np.array([1.0, 0.0, 0.0])

    # ------------------------------------------------------------
    # 1D angular scattering profile (UNCHANGED)
    # ------------------------------------------------------------
    theta_deg, counts_1D, bins, theta_signed_all = angular_scattering_profile(
        vel, norm, tangent=tang, nbins=n_theta, degrees=True
    )



    theta_centers_rad = np.radians(theta_deg)

    # ------------------------------------------------------------
    # φ definition
    # ------------------------------------------------------------
    phi = np.arctan2(vz, vx)
    phi = np.mod(phi, 2*np.pi)

    phi_bins = np.linspace(0, 2*np.pi, n_phi + 1)

    # Match original outputs
    theta_bins = np.linspace(-np.pi/2, np.pi/2, n_theta + 1)
    theta_mid  = 0.5 * (theta_bins[:-1] + theta_bins[1:])

    # Allocate outputs (same shape as before)
    C2D = np.zeros((n_phi, n_theta))
    E2D = np.zeros((n_phi, n_theta))

    # ------------------------------------------------------------
    # Loop over φ bins → call angular_scattering_profile
    # ------------------------------------------------------------
    for j in range(n_phi):
        mask = (phi >= phi_bins[j]) & (phi < phi_bins[j+1])

        if not np.any(mask):
            continue

        vel_subset = vel[mask]

        # --- reuse your 1D function ---
        theta_deg_local, counts_local, bins, theta_signed = angular_scattering_profile(
            vel_subset, norm, tangent=tang, nbins=n_theta, degrees=True
        )



        C2D[j, :] = counts_local

        # --------------------------------------------------------
        # Energy per θ bin (consistent with same θ definition)
        # --------------------------------------------------------
        vx_s, vy_s, vz_s = vel_subset[:, 0], vel_subset[:, 1], vel_subset[:, 2]
        E_particle = 0.5 * (vx_s**2 + vy_s**2 + vz_s**2)

        theta_idx = np.digitize(theta_signed, theta_bins) - 1
        theta_idx = np.clip(theta_idx, 0, n_theta - 1)

        # accumulate (vectorized)
        E_accum = np.zeros(n_theta)
        counts_energy = np.zeros(n_theta)

        np.add.at(E_accum, theta_idx, E_particle)
        np.add.at(counts_energy, theta_idx, 1)

        mask_local = counts_energy > 0
        E_accum[mask_local] /= counts_energy[mask_local]

        E2D[j, :] = E_accum

    return counts_1D, theta_centers_rad, C2D, E2D, theta_mid, phi_bins

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

def compute_drag(datapath):
    import numpy as np
    # Load data
    data = np.loadtxt(datapath, skiprows=9)
    # Extract coordinates and rhon (particle count in cell)
    fx = data[:, 12]
    fy = data[:, 13]
    fz = data[:, 14]
    v1x = data[:, 1]
    v1y = data[:, 2]
    v1z = data[:, 3]
    v2x = data[:, 4]
    v2y = data[:, 5]
    v2z = data[:, 6]
    v3x = data[:, 7]
    v3y = data[:, 8]
    v3z = data[:, 9]

    # px = data[:, 15]
    # py = data[:, 16]
    # pz = data[:, 17]
    # shx = data[:, 18]
    # shy = data[:, 19]
    # shz = data[:, 20]
    area = data[:,23]

    # fx2 = (px + shx)*area

    # exclude vertical sides from total forces computation

    m1 = (v1x == v2x) & (v1x == v3x)                    # vertical, x normal
    m2 = (v1z == v2z) & (v1z == v3z)                    # vertical, z normal
    mv = m1 | m2
    m3 = (v1x == np.min(v1x)) | (v1x == np.max(v1x))    # x edges
    m4 = (v1z == np.min(v1z)) | (v1z == np.max(v1z))    # z edges
    me = m3 | m4
    m = mv & me

    v1 = np.stack([v1x[~m], v1y[~m], v1z[~m]], axis=-1)
    v2 = np.stack([v2x[~m], v2y[~m], v2z[~m]], axis=-1)
    v3 = np.stack([v3x[~m], v3y[~m], v3z[~m]], axis=-1)

    normal = [math.sin(math.radians(70)), math.cos(math.radians(70)), 0]
    # a_norm, angle = projected_triangle_areas(v1, v2, v3, normal)
    # An = np.sum(a_norm, axis=0)
    An = 0
    angle = 0

    ang_x, ang_y, ang_z, azim_xy, azim_yz, azim_zx = face_normal_angles(v1, v2, v3)
    nn = ang_x.shape[0]
    file = datapath.replace("/prop_surfos.1000.dat", "")
    with open("angles"+file+".txt", "w") as writefile:
        for j in range(nn):
            line = f"{ang_x[j]}\t" + f"{ang_y[j]}\t" + f"{ang_z[j]}\t" + f"{azim_xy[j]}\t" + f"{azim_yz[j]}\t" + f"{azim_zx[j]}\t" + "\n"
            writefile.write(line)

    # Output excluded edges to test
    # write_ascii_stl('vertical_faces.stl', v1x[m], v1y[m], v1z[m], v2x[m], v2y[m], v2z[m], v3x[m], v3y[m], v3z[m], solid_name="mesh")

    Fx = np.sum(fx[~m], axis=0)
    Fy = np.sum(fy[~m], axis=0)
    Fz = np.sum(fz[~m], axis=0)
    A = np.sum(area[~m], axis=0)

    return Fx, Fy, Fz, A, An, np.mean(angle)

def write_ascii_stl(filename, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, solid_name="mesh"):
    v1 = np.stack([v1x, v1y, v1z], axis=-1)
    v2 = np.stack([v2x, v2y, v2z], axis=-1)
    v3 = np.stack([v3x, v3y, v3z], axis=-1)

    def compute_normal(a, b, c):
        n = np.cross(b - a, c - a)
        n_norm = np.linalg.norm(n, axis=-1, keepdims=True)
        return n / np.where(n_norm == 0, 1, n_norm)  # Avoid divide-by-zero

    normals = compute_normal(v1, v2, v3)

    with open(filename, 'w') as f:
        f.write(f"solid {solid_name}\n")
        for i in range(len(v1)):
            n = normals[i]
            f.write(f"  facet normal {n[0]} {n[1]} {n[2]}\n")
            f.write("    outer loop\n")
            f.write(f"      vertex {v1[i][0]} {v1[i][1]} {v1[i][2]}\n")
            f.write(f"      vertex {v2[i][0]} {v2[i][1]} {v2[i][2]}\n")
            f.write(f"      vertex {v3[i][0]} {v3[i][1]} {v3[i][2]}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")



def face_normal_angles(v1, v2, v3):
    """
    Compute direction angles of face normals and their azimuths in all 3 planes.

    Parameters
    ----------
    v1, v2, v3 : ndarray (N, 3)
        Vertices of the faces.

    Returns
    -------
    ang_x, ang_y, ang_z : ndarray (N,)
        Angles (in radians) between the normal and the x, y, z axes.
    azim_xy, azim_yz, azim_zx : ndarray (N,)
        Azimuth angles (in radians) in the XY, YZ, and ZX planes:
        - azim_xy: angle from +Y axis to projection in XY plane
        - azim_yz: angle from +Z axis to projection in YZ plane
        - azim_zx: angle from +X axis to projection in ZX plane
    """
    v1, v2, v3 = map(np.asarray, (v1, v2, v3))

    # Compute normals
    e1 = v2 - v1
    e2 = v3 - v1
    n = np.cross(e1, e2)

    # Normalize
    norms = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / norms

    # Clip for numerical safety
    n = np.clip(n, -1.0, 1.0)

    # Direction angles
    ang_x = np.arccos(n[:, 0])
    ang_y = np.arccos(n[:, 1])
    ang_z = np.arccos(n[:, 2])

    # Azimuths
    azim_xy = np.arctan2(n[:, 0], n[:, 1])  # in xy-plane, from +y
    azim_yz = np.arctan2(n[:, 1], n[:, 2])  # in yz-plane, from +z
    azim_zx = np.arctan2(n[:, 2], n[:, 0])  # in zx-plane, from +x

    return ang_x, ang_y, ang_z, azim_xy, azim_yz, azim_zx


def projected_triangle_areas(v1, v2, v3, projection_normal):
    # Normalize projection normal
    n_proj = np.asarray(projection_normal, dtype=np.float64)
    n_proj /= np.linalg.norm(n_proj)

    # Triangle edge vectors
    e1 = v2 - v1
    e2 = v3 - v1

    # Area vectors (triangle normals × 0.5)
    area_vectors = 0.5 * np.cross(e1, e2)
    
    # Projected area: |dot(n_proj, area_vector)|
    projected_areas = np.abs(np.dot(area_vectors, n_proj))

    # Normalize triangle normals
    norms = np.linalg.norm(area_vectors, axis=1, keepdims=True)
    normals = area_vectors / np.where(norms == 0, 1, norms)

    # Project both triangle normals and projection normal into XY plane
    normals_xy = normals[:, :2]  # drop Z
    proj_xy = n_proj[:2]         # drop Z

    # Normalize projected vectors
    norms_xy = np.linalg.norm(normals_xy, axis=1, keepdims=True)
    proj_xy_norm = np.linalg.norm(proj_xy)

    # Prevent division by zero
    proj_xy_safe = proj_xy / (proj_xy_norm if proj_xy_norm != 0 else 1)
    normals_xy_safe = normals_xy / np.where(norms_xy == 0, 1, norms_xy)

    # Compute dot product between triangle normals and projection vector in XY
    dot_xy = np.clip(np.sum(normals_xy_safe * proj_xy_safe, axis=1), -1.0, 1.0)
    
    # Angle in degrees
    angles_xy_deg = np.degrees(np.arccos(dot_xy))

    return projected_areas, angles_xy_deg


def vector_angle(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Compute the angle between rows of two N×3 arrays in degrees.

    Parameters
    ----------
    A : np.ndarray
        N×3 array where each row is a 3D vector.
    B : np.ndarray
        N×3 array where each row is a 3D vector.

    Returns
    -------
    theta_deg : np.ndarray
        N×1 array of angles in degrees.
    """
    # Validate shape
    if A.shape != B.shape:
        raise ValueError("Input arrays must have the same shape.")

    # Compute dot products for each row
    dotAB = np.sum(A * B, axis=1)

    # Compute norms of each vector
    normA = np.linalg.norm(A, axis=1)
    normB = np.linalg.norm(B, axis=1)

    # Compute cosine of angle, clamp to [-1, 1] for numerical safety
    cos_theta = dotAB / (normA * normB)
    # cos_theta = np.clip(cos_theta, -1.0, 1.0)

    # Compute angle in degrees
    theta_deg = np.degrees(np.arccos(cos_theta))

    # Assign sign based on x-component of A
    theta_deg = theta_deg * np.sign(A[:, 0])

    return theta_deg


    
def vector_angle_projected(
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray
) -> np.ndarray:
    """
    Signed angle between A and the projection of B onto the plane
    defined by A and C (right-hand rule about n = A × C).

    Parameters
    ----------
    A, B, C : np.ndarray
        N×3 arrays of vectors.

    Returns
    -------
    theta_deg : np.ndarray
        Signed angles in degrees, shape (N,)
    """
    if A.shape != B.shape or A.shape != C.shape:
        raise ValueError("Input arrays must all have the same shape.")

    # Plane normal
    n = np.cross(A, C)
    n_norm = np.linalg.norm(n, axis=1, keepdims=True)
    if np.any(n_norm == 0.0):
        raise ValueError("A and C must not be collinear.")

    n_hat = n / n_norm

    # Project B onto plane
    B_proj = B - np.sum(B * n_hat, axis=1, keepdims=True) * n_hat

    # Normalize A and projected B
    A_hat = A / np.linalg.norm(A, axis=1, keepdims=True)
    Bp_hat = B_proj / np.linalg.norm(B_proj, axis=1, keepdims=True)

    # Unsigned angle
    cos_theta = np.sum(A_hat * Bp_hat, axis=1)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)

    # Sign from right-hand rule about plane normal
    cross_ABp = np.cross(A_hat, Bp_hat)
    sign = np.sign(np.sum(cross_ABp * n_hat, axis=1))

    theta_signed_deg = np.degrees(theta * sign)
    return theta_signed_deg




main()