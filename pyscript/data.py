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

ROOT = f'../bin/'
tmax        = 0.40          # total time    
dtoutput    = 0.01          # dt for output 
nstep = int(tmax/dtoutput)
timelist = np.linspace(0, tmax, nstep+1)

plt.figure(figsize=(8,6))

for frame in range(0,nstep+1,8):
    filename = ROOT + f"output_"+str(frame).zfill(5)
    pts = read_particles(filename)
    select = np.sqrt(pts[:,X2]**2+pts[:,X3]**2)<=0.4
    r = np.sqrt(pts[select,X1]**2 + pts[select,X2]**2 + pts[select,X3]**2)
    plt.scatter(r, pts[select,DEN], s = .05)
plt.xlim(0, 2.5)
plt.xlabel("r", fontsize = 15)
plt.ylabel(r'$\rho$', fontsize = 15)
plt.savefig("density-r.png", bbox_inches='tight', dpi=300)
# plt.show()
plt.close()

def cal_total_e(pts):
    ie = 0 
    ke = 0
    tote = 0
    for i in range(len(pts[:,X1])):
        ie = ie + pts[i,MAS]*pts[i, ENE]
        ke = ke + 0.5*pts[i,MAS]*(pts[i,VEL1]**2 + pts[i,VEL2]**2 + pts[i,VEL3]**2)
        # print(ke)
    tote = ie+ke 
    return ie, ke, tote 

ie = np.zeros(len(timelist))
ke = np.zeros(len(timelist))
tote = np.zeros(len(timelist))
for frame in range(0, nstep+1):
    filename = ROOT + f"output_"+str(frame).zfill(5)
    pts = read_particles(filename)
    ie[frame], ke[frame], tote[frame] = cal_total_e(pts) 

plt.figure(figsize=(8,6))
plt.plot(timelist, tote, label = 'total e')
plt.plot(timelist, ie, label = 'ie')
plt.plot(timelist, ke, label = 'ke')
plt.legend(fontsize = 15)
plt.xlabel('t', fontsize = 15)
plt.ylabel('E', fontsize = 15)
# plt.show()
plt.savefig("energy-t.png", bbox_inches='tight', dpi=300)
plt.close()