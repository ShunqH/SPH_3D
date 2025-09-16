import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt 
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

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

def kernel(r, h):
    r = np.atleast_1d(r)
    h = np.atleast_1d(h)
    q = r / h
    W = np.zeros_like(q)
    factor = 1.0 / (np.pi * h**3)

    mask1 = q < 1.0
    mask2 = (q >= 1.0) & (q < 2.0)

    W[mask1] = factor[mask1] * (1 - 1.5*q[mask1]**2 + 0.75*q[mask1]**3)
    W[mask2] = factor[mask2] * 0.25*(2 - q[mask2])**3

    return W if len(W) > 1 else W.item()


def render_den(x, y, z, pts, hfact = 3):
    nx, ny, nz = len(x), len(y), len(z)
    floor = 1e-8
    rho = np.ones((nz, ny, nx)) * floor
    
    gas = pts[:, STA] >= -10
    coords = np.column_stack((pts[gas, X1], pts[gas, X2], pts[gas, X3]))
    mass = pts[gas, MAS]
    h    = pts[gas, LEN]

    tree = cKDTree(coords)
    dr = np.max(h)

    XX, YY, ZZ = np.meshgrid(x, y, z, indexing="ij")
    grid = np.column_stack((XX.ravel(), YY.ravel(), ZZ.ravel()))

    neighbors_list = tree.query_ball_point(grid, r=hfact*dr)

    rho_flat = np.zeros(len(grid)) + floor
    for gi, neighbors in enumerate(neighbors_list):
        if not neighbors:
            continue
        r = grid[gi]
        neigh_coords = coords[neighbors]
        neigh_h = h[neighbors]
        neigh_m = mass[neighbors]

        dist = np.linalg.norm(neigh_coords - r, axis=1)
        mask = dist <= hfact*neigh_h
        if np.any(mask):
            W = kernel(dist[mask], neigh_h[mask])
            rho_flat[gi] += np.sum(neigh_m[mask] * W)

    rho = rho_flat.reshape(nx, ny, nz).transpose(2,1,0)
    return rho

def GridMesh(xmin, xmax, n):
    dx = (xmax-xmin)/n 
    return np.linspace(xmin+dx/2, xmax-dx/2, n)


ROOT = f'../bin/sod/'
tmax        = 0.50          # total time    
dtoutput    = 0.01          # dt for output 
nstep = int(tmax/dtoutput)
timelist = np.linspace(0, tmax, nstep+1)

# plot 1D density profile
plt.figure(figsize=(8,6))
for frame in range(0,nstep+1,10):
    filename = ROOT + f"output_"+str(frame).zfill(5)
    pts = read_particles(filename)
    select = np.sqrt(pts[:,X2]**2+pts[:,X3]**2)<=10
    order = np.argsort(pts[select,X1])
    plt.scatter(pts[select,X1], pts[select,DEN], s = 0.1)
    plt.plot(pts[select,X1][order], pts[select,DEN][order], 
             label=f"t={frame*dtoutput:.2f}")
plt.xlim(-1, 1)
plt.xlabel("x", fontsize = 15)
plt.ylabel(r'$\rho$', fontsize = 15)
plt.savefig("density-r.png", bbox_inches='tight', dpi=300)
# plt.show()
plt.close()

# plot energy time evolution
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
plt.hlines(y=tote[0],xmin=timelist[0],xmax=timelist[-1],color='gray')
plt.legend(fontsize = 15)
plt.xlabel('t', fontsize = 15)
plt.ylabel('E', fontsize = 15)
plt.xlim(timelist[0], timelist[-1])
plt.tick_params(axis='both', labelsize=15, which='both', top=True, bottom=True, left=True, right=True)  
# plt.show()
plt.savefig("energy-t.png", bbox_inches='tight', dpi=300)
plt.close()

def add_colorbar(fig, ax, pcm):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.0)
    cbar = fig.colorbar(pcm, cax=cax, extend='neither')
    # cbar.ax.set_yscale('log')
    return cbar

# render density in mesh grid, 
# set nz = 1 to render cross-section at z=0
for frame in range(0, 51, 50):
    filename = ROOT + f"output_"+str(frame).zfill(5)
    pts = read_particles(filename)
    x = GridMesh(-1.0, 1.0, 100)
    y = GridMesh(-0.5, 0.5, 50)
    z = GridMesh(-0.1, 0.1, 1)
    z_index = 0
    rho = render_den(x, y, z, pts)
    XX, YY = np.meshgrid(x,y)
    fig, ax = plt.subplots(1, 1, sharey=True, figsize=(8,6))
    pcm = plt.pcolormesh(XX, YY, rho[z_index,:,:], 
                norm=colors.Normalize(vmin=0, vmax=1),
                cmap = "gnuplot"
                )
    plt.gca().set_aspect('equal')
    cbar = add_colorbar(fig, ax, pcm)
    plt.savefig(f"Mesh_rho_"+str(frame).zfill(5)+f".png", bbox_inches='tight', dpi=300)
    plt.close() 
    # plt.show()