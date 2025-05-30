import numpy as np
import matplotlib.pyplot as plt 

X1 = 0      # position x1
VEL1 = 1    # velocity in x1 direction 
MAS = 2     # mass 
LEN = 3     # smoothing length
DEN = 4     # density
ENE = 5     # internal energy (u) 
PRE = 6     # pressure 
SPE = 7     # sound speed
ACC1 = 8    # acceleration for velosity
ACCU = 9    # acceleration for internal energy

def read_particles(filename):
    with open(filename, 'rb') as f:
        rows = np.fromfile(f, dtype=np.int32, count=1)[0]  
        cols = np.fromfile(f, dtype=np.int32, count=1)[0] 

        particles = np.fromfile(f, dtype=np.float64) 

        particles = particles.reshape((rows, cols))
    
    return particles

plt.figure(figsize = (12,8))
path = "../bin/"
dt = 0.01
nghost = 6

for i in range(0, 21, 4):
    filename = path + 'output_' + str(i).zfill(5)
    particles = read_particles(filename)
    t = i*dt 
    act_par = particles[nghost:-nghost]
    #plt.scatter(act_par[:, X1], act_par[:, DEN], label = "t = "+ format(t, ".3f"), s = 5)
    plt.plot(act_par[:, X1], act_par[:, DEN], label = "t = "+ format(t, ".3f"), marker = "o")
#plt.xlim(-0.05,0.05)
plt.grid()
plt.legend(fontsize = 18)
plt.xlabel("x", fontsize = 18)
plt.ylabel(r"$\rho$", fontsize = 18)
plt.savefig("sod.png", bbox_inches='tight', dpi=300)
plt.close()