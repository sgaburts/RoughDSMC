import shutil
import numpy as np

filename = "new_template.in"
change_lines = np.array([14, 22, 24, 38, 40, 84])
change_lines = change_lines-1

theta_beam = np.linspace(0, 90, 10)
# xy_calibration = np.array([5.663, 2.693, 1.368, 0.686, 0.276, 0.095])
xy_calibration = 2.693
# fnum = np.array([1e10, 5e8, 3e7, 2e6, 5e4, 7e2])
fnum = "1e11" #"1e-2"
# surf_name = "Ge_Si_Exposed_Height_Processed"
# surf_name = "Ge_Si_Height_Processed"
surf_name = "smooth_plate_thin"

y_scale = 1

# Not used
scale = 1e-3
nx = 1024
nz = 768

max_y_scatter = 0.2 #0.0001  # max height to dump particles for scattering

n = np.size(theta_beam)

velocity = 8000
temp = 0.0001

scattering = 'diffuse' # diffuse, specular or cll

surf_temp = 300 # for diffuse or cll

acc_nor = 0.30  # normal for CLL
acc_tan = 0.58  # tangential for CLL


for i in range(n):
    new_file = scattering+"_"+str(int(theta_beam[i]))+"deg.in"
    # new_file = scattering+"_"+str(y_scale)+"ys_"+str(int(theta_beam[i]))+"deg.in"

    # Not used
    xp = 0.5*xy_calibration*scale*nx
    zp = 0.5*xy_calibration*scale*nz

    xlh = 1.18*xp
    yl  = 1.18*xp
    zlh = zp + 0.01*zp
    

    new_line_content = [""] * 6 
    # new_line_content[0] = "create_box          -"+str(xlh)+" "+str(xlh)+" 0 "+str(yl)+" -"+str(zlh)+" "+str(zlh)+"\n"
    # new_line_content[0] = "create_box          -0.0002 0.0002 0 0.0002 -0.000133 0.000133\n"
    new_line_content[0] = "create_box          -0.55 0.55 0 0.6 -0.55 0.55\n"
    new_line_content[1] = "global              nrho 1.46e15 fnum "+fnum+"\n"  
    new_line_content[2] = "mixture             atm N2 O2 O vstream "+str(velocity*np.sin(np.radians(theta_beam[i])))+" -"+str(velocity*np.cos(np.radians(theta_beam[i])))+" 0.0 temp "+str(temp)+" # "+str(int(theta_beam[i]))+"deg\n"

    # new_line_content[3] = "read_surf           "+surf_name+" group rocket trans 0 0 0\n"  
    new_line_content[3] = "read_surf           "+surf_name+" group rocket trans 0 0 0 scale 1 "+str(y_scale)+" 1\n"  


    if scattering == 'specular':
        new_line_content[4] = "surf_collide        1 specular # new\n"
    if scattering == 'diffuse':
        new_line_content[4] = "surf_collide        1 cll "+str(surf_temp)+" 1.0 1.0 0.0 0.0\n"
    if scattering == 'cll':
        new_line_content[4] = "surf_collide        1 cll "+str(surf_temp)+" "+str(acc_nor)+" "+str(acc_tan)+" 0.0 0.0\n"

    new_line_content[5] = "region             scatter_box block -INF INF  0.0 "+str(max_y_scatter)+"  -INF INF\n"  
    


    # Step 1: Make a backup
    shutil.copyfile(filename, new_file)

    # Step 1: Read all lines
    with open(new_file, 'r') as file:
        lines = file.readlines()

    # Step 2: Modify the specific line  
    for j in range(np.size(change_lines)):
        if 0 <= change_lines[j] < len(lines):
            lines[change_lines[j]] = new_line_content[j]
        else:
            print("Line number out of range.")

    # Step 3: Write lines back to the file
    with open(new_file, 'w') as file:
        file.writelines(lines)
