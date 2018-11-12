import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm

## Rough first pass. Assume form of grid is 13 rows:
## X  Y  Z  DEN ABUN TEM DTEM V_TURB V_X V_Y V_Z H_FRAC H2_FRAC

instr = raw_input('Gridfile to visualize?')

with open(instr,'r') as f:
	grid = np.loadtxt(f)

z = grid.transpose()[2]/1.496e11

r_sph = ((grid.transpose()[0]/1.496e11)**2+
        (grid.transpose()[1]/1.496e11)**2+
         (z)**2)**0.5

r_pol = ((grid.transpose()[0]/1.496e11)**2+
        (grid.transpose()[1]/1.496e11)**2)**0.5

print np.min(z),np.max(z)

print np.min(r_sph),np.max(r_sph)
print np.min(r_pol),np.max(r_pol)

print instr.split('.')[0]
for colorcol in range(3,7)+range(11,13):
    print colorcol
    color = 1.#np.log10(grid[:,colorcol])
    plt.figure(figsize=(20,10))
#    plt.scatter(grid.transpose()[0]/1.496e11,grid.transpose()[2]/grid.transpose()[0],
#    plt.scatter(grid.transpose()[0]/1.496e11,grid.transpose()[2]/1.496e11,
    plt.scatter(r_pol,z/r_pol,
    			c=np.log10(grid.transpose()[colorcol]),marker='.')
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim([.08,2])
    #plt.ylim([0.01,0.8])
    plt.colorbar()
    #plt.savefig(instr.split('.')[0]+'_col'+str(colorcol)+'.png',fmt='png',dpi=300)