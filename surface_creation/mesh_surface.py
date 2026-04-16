import numpy as np
import trimesh
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

def main():

    nx = 100 #512
    nz = 100 #512

    lx = 1 #nx*0.095*1e-3
    lz = 1 #nz*0.095*1e-3
    ly = 0.000001

    targ = 0.001

    scale = 1 #1e-3

    subsample = 8 #8 #4 #subsampling for random rms

    filename1 = "smooth_plate.stl"
    filename2 = "test"
    surfacedata = "Ge_Kapton_30nm_After_exposure_Height.csv"

    x = np.linspace(-0.5*lx, 0.5*lx, nx)
    z = np.linspace(-0.5*lz, 0.5*lz, nz)

    generate_plate_stl(lx=lx,ly=ly,lz=lz,nx=nx,nz=nz,filename=filename1)
    print('Smooth plate generated')
    
    # Load the mesh
    mesh = trimesh.load_mesh(filename1)

    # Load heightmap
    # height_map = read_custom_csv(surfacedata)

    # Random with given RMS
    height_map = generate_random_matrix_with_rms((nx//subsample, nz//subsample), target_rms=targ)
    # print("RMS:", np.sqrt(np.mean(height_map**2)))

    # Scale
    # height_map = height_map*scale + 1e-3
    height_map = height_map*scale
    height_map = height_map + abs(np.min(height_map))
    print(height_map.shape)

    final_rms = np.sqrt(np.sum((height_map-abs(np.min(height_map)))**2)/(nx*nz))
    # print("RMS:", final_rms)
    # print("RMS:", np.sqrt(np.mean(height_map**2)))

    # Example height map
    # nx, nz = 60, 40
    # x = np.linspace(-0.005, 0.005, nx)
    # z = np.linspace(-0.0025, 0.0025, nz)
    # X, Z = np.meshgrid(x, z)
    # height_map = np.random.uniform(0.0001, 0.0002, size=(nx, nz))
    # nx, nz = 1024, 768
    # x = np.linspace(-0.5*nx, 0.5*nx, nx)
    # z = np.linspace(-0.5*nz, 0.5*nz, nz)
    # X, Z = np.meshgrid(x, z)
    # height_map = np.random.uniform(1, 6, size=(nx, nz))

    # Plot the height map
    plt.figure(figsize=(6, 5))
    plt.imshow(height_map.T,
            origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar(label='Height (y)')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title('Height Map Applied to Top Surface (in y)')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(f"heightmap.png", dpi=300)

    # Ranges for the interpolator
    x_range = (x.min(), x.max())
    z_range = (z.min(), z.max())

    # Apply height map
    mesh = update_top_surface_y(mesh, height_map, x_range, z_range)

    print('Plate deformed')

    # Cleanup
    mesh.remove_unreferenced_vertices()
    mesh.update_faces(mesh.unique_faces())
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.fix_normals()

    print('Mesh cleaned')

    # Export updated mesh
    # mesh.export(file_obj=filename2+str(targ)+'.stl', file_type='stl_ascii')
    mesh.export(file_obj=filename2+'.stl', file_type='stl_ascii')

    

#=============================================================================================

def read_custom_csv(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find the line that contains just "Height"
    for i, line in enumerate(lines):
        if line.strip().strip('"') == "Height":
            data_start = i + 1
            break
    else:
        raise ValueError("No 'Height' label found.")

    # Read the rest of the lines as data
    data_lines = lines[data_start:]
    data = []
    for line in data_lines:
        row = [float(val.strip().strip('"')) for val in line.strip().split(',')]
        data.append(row)

    return np.array(data)

def generate_plate_stl(lx, ly, lz, nx, nz, filename):
    """
    Generate a watertight plate mesh as an ASCII STL file.

    Plate spans x ∈ [-0.5*lx, 0.5*lx], y ∈ [0, ly], z ∈ [-0.5*lz, 0.5*lz]
    with resolution nx in x, nz in z, and 2 layers in y (top and bottom).

    Parameters:
        lx, ly, lz: Plate dimensions.
        nx, nz: Resolution in x and z directions.
        filename: Output filename.
    """
    ny = 2
    x = np.linspace(-0.5 * lx, 0.5 * lx, nx)
    y = np.array([0, ly])
    z = np.linspace(-0.5 * lz, 0.5 * lz, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    vertices = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    # vertices = np.round(vertices, decimals=6)

    def idx(i, j, k): return i * ny * nz + j * nz + k

    faces = []

    # --- Top and bottom surfaces ---
    for i in range(nx - 1):
        for k in range(nz - 1):
            # Bottom (y = 0): j = 0
            v00 = idx(i, 0, k)
            v10 = idx(i + 1, 0, k)
            v11 = idx(i + 1, 0, k + 1)
            v01 = idx(i, 0, k + 1)
            faces += [[v00, v10, v11], [v00, v11, v01]]

            # Top (y = ly): j = 1
            v00 = idx(i, 1, k)
            v10 = idx(i + 1, 1, k)
            v11 = idx(i + 1, 1, k + 1)
            v01 = idx(i, 1, k + 1)
            faces += [[v00, v11, v10], [v00, v01, v11]]  # flipped for outward normal

    # --- Side walls ---
    for j in [0, 1]:  # bottom and top layers
        for i in range(nx - 1):
            # Front face (z = 0)
            v0 = idx(i, j, 0)
            v1 = idx(i + 1, j, 0)
            v2 = idx(i + 1, 1 - j, 0)
            v3 = idx(i, 1 - j, 0)
            if j == 0:  # bottom layer
                faces += [[v0, v1, v2], [v0, v2, v3]]
        for k in range(nz - 1):
            # Right face (x = +)
            v0 = idx(nx - 1, j, k)
            v1 = idx(nx - 1, j, k + 1)
            v2 = idx(nx - 1, 1 - j, k + 1)
            v3 = idx(nx - 1, 1 - j, k)
            if j == 0:
                faces += [[v0, v1, v2], [v0, v2, v3]]
        for i in range(nx - 1):
            # Back face (z = +)
            v0 = idx(i, j, nz - 1)
            v1 = idx(i + 1, j, nz - 1)
            v2 = idx(i + 1, 1 - j, nz - 1)
            v3 = idx(i, 1 - j, nz - 1)
            if j == 0:
                faces += [[v0, v2, v1], [v0, v3, v2]]
        for k in range(nz - 1):
            # Left face (x = -)
            v0 = idx(0, j, k)
            v1 = idx(0, j, k + 1)
            v2 = idx(0, 1 - j, k + 1)
            v3 = idx(0, 1 - j, k)
            if j == 0:
                faces += [[v0, v2, v1], [v0, v3, v2]]

    mesh = trimesh.Trimesh(vertices=vertices, faces=np.array(faces), process=False)
    # Clean mesh (future-proof replacements)
    mesh.update_faces(mesh.unique_faces())
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    mesh.fix_normals()
    # To ensure watertightness for SPARTA
    mesh = mesh.process(validate=True)

    # Export as ASCII STL
    mesh.export(file_obj=filename, file_type='stl_ascii')
    print(f"Watertight STL saved to: {filename}")

def update_top_surface_y(mesh, height_map, x_range, z_range):
    """
    Updates the y-coordinate (height) of the top surface using a height map defined over (x, z).

    Parameters:
        mesh (trimesh.Trimesh): The loaded STL mesh.
        height_map (np.ndarray): 2D array of new y-values (heights) with shape (nz, nx).
        x_range (tuple): (min_x, max_x) corresponding to axis-1 of height_map.
        z_range (tuple): (min_z, max_z) corresponding to axis-0 of height_map.
    """
    vertices = mesh.vertices.copy()

    # Get shape and create axis vectors
    nz, nx = height_map.shape
    x_vals = np.linspace(x_range[0], x_range[1], nx)
    z_vals = np.linspace(z_range[0], z_range[1], nz)

    # Interpolator over (z, x)
    interpolator = RegularGridInterpolator((z_vals, x_vals), height_map, bounds_error=False, fill_value=None)

    # Threshold to detect the top surface
    y_threshold = np.percentile(vertices[:, 1], 90)

    for i, (x, y, z) in enumerate(vertices):
        if y >= y_threshold:
            interpolated_y = interpolator((z, x))  # Notice: (z, x)
            if interpolated_y is not None:
                vertices[i][1] = interpolated_y

    mesh.vertices = vertices
    return mesh

def generate_random_matrix_with_rms(shape, target_rms, distribution='normal'):
    if distribution == 'normal':
        A = np.random.randn(*shape)  # standard normal distribution
    elif distribution == 'uniform':
        A = np.random.rand(*shape) - 0.5  # uniform in [-0.5, 0.5]
    else:
        raise ValueError("Unsupported distribution")

    # current_rms = np.sqrt(np.mean(A**2))
    current_rms = np.sqrt(np.sum(A**2)/(shape[0]*shape[1]))
    scaling_factor = target_rms / current_rms
    A_scaled = A * scaling_factor

    final_rms = np.sqrt(np.sum(A_scaled**2)/(shape[0]*shape[1]))
    print("RMS:", final_rms)

    return A_scaled

main()



