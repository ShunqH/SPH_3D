import numpy as np
import matplotlib.pyplot as plt 

X1      = 0     # position x
X2      = 1     # position y
X3      = 2     # position z
VEL1    = 3     # velocity in x direction 
VEL2    = 4     # velocity in y direction 
VEL3    = 5     # velocity in z direction 
ACC1    = 6     # velocity in x direction 
ACC2    = 7     # velocity in y direction 
ACC3    = 8     # velocity in z direction 
ACCU    = 9     # velocity in z direction 
DEN     = 10    # density
MAS     = 11    # mass 
LEN     = 12    # smoothing length
ENE     = 13    # internal energy (u) 
STA     = 14    # status. 0: alive; 1: dead 

def read_particles(filename):
    with open(filename, 'rb') as f:
        n = np.fromfile(f, dtype=np.int32, count=1)[0]
        x1    = np.fromfile(f, dtype=np.float64, count=n)
        x2    = np.fromfile(f, dtype=np.float64, count=n)
        x3    = np.fromfile(f, dtype=np.float64, count=n)
        vel1  = np.fromfile(f, dtype=np.float64, count=n)
        vel2  = np.fromfile(f, dtype=np.float64, count=n)
        vel3  = np.fromfile(f, dtype=np.float64, count=n)
        acc1  = np.fromfile(f, dtype=np.float64, count=n)
        acc2  = np.fromfile(f, dtype=np.float64, count=n)
        acc3  = np.fromfile(f, dtype=np.float64, count=n)
        accu  = np.fromfile(f, dtype=np.float64, count=n)
        den   = np.fromfile(f, dtype=np.float64, count=n)
        mas   = np.fromfile(f, dtype=np.float64, count=n)
        length= np.fromfile(f, dtype=np.float64, count=n)
        ene   = np.fromfile(f, dtype=np.float64, count=n)
        status= np.fromfile(f, dtype=np.int32,   count=n)
    gamma = 1.4
    pres = (gamma-1.)*den*ene
    data = np.stack([x1, x2, x3, vel1, vel2, vel3, acc1, acc2, acc3, accu, den, mas, length, ene, status], axis=1)

    return data

pts = read_particles(f"../bin/output_00000")
select = np.sqrt(pts[:,X2]**2+pts[:,X3]**2)<=0.2
plt.figure(figsize=(8,6))
for frame in range(0,21,4):
    filename = f"../bin/output_"+str(frame).zfill(5)
    pts = read_particles(filename)
    plt.scatter(pts[select,X1], pts[select,DEN], s = .05)
plt.xlim(-1.2, 1.2)
plt.savefig("figure.png", bbox_inches='tight', dpi=300)
plt.show()