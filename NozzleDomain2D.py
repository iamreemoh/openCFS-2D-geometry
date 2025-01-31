import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from mesh_tool import *  
from scipy.interpolate import interp1d

R_I = 1.0       # Inlet radius
R_star = 0.5    # Throat radius
r_star = 0.20
theta = np.radians(15)  # Divergence angle in radians
L = 5.0       # Total nozzle length

d = (3 / 2) * R_star * np.tan(theta)
x_plus = L - d
r_plus = (5 * d**2) / (12 * R_star) + r_star
a = 0.2  # Cylindrical section length 

b = 2 * (x_plus + ((r_plus - R_I) / np.tan(theta)) - a)
e = (R_I - r_star) / np.tan(theta)
c = e - (b / 2) - (5 * d / 8)

x_values = np.linspace(0, L, 1000)  # Divide length into 1000 segments
r_values = np.zeros_like(x_values)  # Initialize radius array

# Calculate radius for each segment
for i, x in enumerate(x_values):
    if 0 <= x <= a:
        r_values[i] = R_I
    elif a <= x <= a + b:
        r_values[i] = R_I - (b * np.tan(theta) / 2) * (((x - a) / b) ** 3) * (2 - ((x - a) / b))
    elif a + b <= x <= a + b + c:
        r_values[i] = R_I + (a + (b / 2)) * np.tan(theta) - x * np.tan(theta)
    else:
        r_values[i] = ((L - x)**2 / (12*R_star))*(6 - ((L - x) / d) ** 2) + r_star

# Reflect the nozzle shape about the vertical line at x = L
x_mirror = 2 * L - x_values  # Reflect x_values about x = L
r_mirror = r_values  # r_values remain the same

# Stretch the mirrored x-values by a factor of 3
x_mirror_stretched = L + 3 * (x_mirror - L)

# Stretch the mirrored r-values by a factor of 3, and subtract the original max radius value to avoid shifting
r_mirror_stretched = 0.2 + 2 * (r_mirror - np.min(r_mirror))  # Subtract the min value to keep symmetry

N = 200
x_mirror_stretched = x_mirror_stretched[N:]
r_mirror_stretched = r_mirror_stretched[N:]

# Reverse the mirrored arrays
x_mirror_stretched_reversed = x_mirror_stretched[::-1]
r_mirror_stretched_reversed = r_mirror_stretched[::-1]

# Concatenate the original arrays with the reversed mirrored arrays
x_temp = np.concatenate((x_values, x_mirror_stretched_reversed))
r_temp = np.concatenate((r_values, r_mirror_stretched_reversed))

# Resolution based on x-resolution
nx = 200 # t needs to match c * width/nx with c = 1,2,3,..
# Remove duplicate x-values
unique_indices = np.unique(x_temp, return_index=True)[1]
x_temp_unique = x_temp[unique_indices]
r_temp_unique = r_temp[unique_indices]

M = 150
interp_func = interp1d(x_temp_unique, r_temp_unique, kind='cubic')
# Generate new x values with only the desired number of points
x_combined_unfiltered = np.linspace(x_temp.min(), x_temp.max(), nx + M)
r_combined_unfiltered = interp_func(x_combined_unfiltered)

x_combined = x_combined_unfiltered[:-M]
r_combined = r_combined_unfiltered[:-M]
# dimensions of box (m)
width = 30
height = 3

# inner box (m)
wi = x_combined[nx -1] - x_combined[0]
hi = r_combined[nx -1]

print("x_combined",len(x_combined))
print("r_combined",len(r_combined))
print("hi",r_combined[nx-1])
print("wi",x_combined[nx-1] - x_combined[0])
# # Resolution based on x-resolution
# nx = 40  # t needs to match c * width/nx with c = 1,2,3,..
ny = int((height / width) * nx)
mesh = create_2d_mesh(nx, ny, width, height)

# Interface region: solid with thickness t (in m)
t = 2 * width / nx  # needs to match the discretization
print("t",t)
step = 0
# Define regions
for e in mesh.elements:
    x, y = mesh.calc_barycenter(e)
    x = round(x, 2)
    y = round(y, 2)

    # Get the corresponding radius value
    if step < len(r_combined):  # Ensure step is within bounds
        r_with_thickness = r_combined[step] + t
        r_with_thickness = round(r_with_thickness, 2)

        # if y >= r_combined[step] and y <= r_with_thickness:
        if y >= r_combined[step] and y <= r_with_thickness:
            e.region = 'solid'
        elif y < r_combined[step]:
            e.region = 'void'
        elif y >= r_combined[nx -1] + 1.5 and y <= r_combined[nx -1] + t + 1.5: # and x <= 14:
            e.region = 'solid'
    step = (step + 1) % len(r_combined)  # Reset step if it reaches the end

box_width = []    
box_height = []
box_curve = []
casing = []  

step = 0    
for i, n in enumerate(mesh.nodes):
  x, y = n
  if x >= x_combined[0] and x <= x_combined[nx-1]: # only works if resolution hits the nodes - in case of issues round()
    box_width.append(i)
  if y >= r_combined[0] and y <= r_combined[nx-1]:
    box_height.append(i)
  if y >= r_combined[nx -1] + 1.5 and y <= r_combined[nx -1] + t + 1.5:# and x <= 14:
    casing.append(i)
  if step == nx:
    step = 0

  elif step < nx and y >= r_combined[step] and y <= r_combined[step] + t and x >= x_combined[0] and x <= x_combined[nx-1]:
#   elif step < nx and y == r_combined[step] and x >= x_combined[0] and x <= x_combined[nx-1]:
    box_curve.append(i)
  step += 1
  

mesh.bc.append(('box_width',box_width))
mesh.bc.append(('box_height',box_height))
mesh.bc.append(('box_curve',box_curve))
mesh.bc.append(('casing',casing))       

print('size of box_width/box_height/box_curve',len(box_width),len(box_height),len(box_curve))    
    
temp = round(t,2)
    
f = 'box2d-t_' + str(temp) + '-nx_' + str(nx) + '-ny_' + str(ny) + '.mesh'   

write_ansys_mesh(mesh, f)
print('created ', f) 
# hint: cfs -m <mesh> <problem> -g will just write the mesh without simulation/optimization for check in ParaView


import matplotlib.pyplot as plt

# Separate nodes into solid, void, and default (mech) regions
solid_region_x, solid_region_y = [], []
void_region_x, void_region_y = [], []
default_region_x, default_region_y = [], []
# Classify nodes by their regions
count = 0
for step, n in enumerate(mesh.nodes):
    x, y = n
    if count == nx:
        count = 0
    elif count < nx:
        r_at_x = r_combined[count]
        r_with_thickness = r_at_x + t

    # Solid region (outer box with thickness t)
    if count < nx and y >= r_at_x and y <= r_with_thickness:
        solid_region_x.append(x)
        solid_region_y.append(y)
    # Void region (inner box)
    elif y < r_at_x and count < nx:
        void_region_x.append(x)
        void_region_y.append(y)
    # Default (mech region)
    else:
        default_region_x.append(x)
        default_region_y.append(y)

    count = (count + 1) % len(r_combined)  # Reset count properly

# Plot the nodes
plt.figure(figsize=(12, 10))

# Scatter plots for each region
plt.scatter(default_region_x, default_region_y, color='gray', label='Mech Region', alpha=0.5, s=10)
plt.scatter(solid_region_x, solid_region_y, color='red', label='Solid Region', alpha=0.8, s=10)
plt.scatter(void_region_x, void_region_y, color='blue', label='Void Region', alpha=0.8, s=10)

# Create gridlines for the discretized mesh
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Add labels and legend
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Discretized Mesh with Regions')
plt.legend()

# Dynamically set axis limits to include all points
all_x = default_region_x + solid_region_x + void_region_x
all_y = default_region_y + solid_region_y + void_region_y
x_margin = 1  # Add margin to the x-axis
y_margin = 1  # Add margin to the y-axis
plt.xlim(min(all_x) - x_margin, max(all_x) + x_margin)
plt.ylim(min(all_y) - y_margin, max(all_y) + y_margin)

# Show the plot
plt.show()


print("RUN: cfs -m ",f," 2D_distributed_load")
