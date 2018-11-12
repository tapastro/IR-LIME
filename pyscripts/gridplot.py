import numpy as np
from matplotlib import pyplot as plt

griddir = '/Users/tpauly/projects_unsync/lime/models/sr21/output/'
#gridname = 'grid_interp_nu2j21_popdiag_flattable.vtk'
gridname = 'grid_sr21_nu2j21_0p9tem.vtk'

with open(griddir+gridname,'r') as f:
	header = [f.readline() for _ in range(4)]
	npts = int(f.readline().split()[1])
	ptarr = np.zeros((npts,3))
	print np.shape(ptarr)
	for i in range(npts):
		ptarr[i] = map(float,f.readline().split())
	blank = f.readline()
	ncells = int(f.readline().split()[1])
	for j in range(ncells):
		f.readline()
	blank = f.readline()
	#celltypes next
	cellheader = f.readline()
	for k in range(ncells):
		f.readline()
	# No separating newline after celltypes?
	#blank = f.readline()
	# H2 Density next
	h2header = [f.readline() for _ in range(3)]

	h2denarr = np.zeros((npts,))
	for m in range(npts):
		h2denarr[m] = float(f.readline())
	# Mol abun next
	molheader = [f.readline() for _ in range(2)]
	molabunarr = np.zeros((npts,))
	for n in range(npts):
		molabunarr[n] = float(f.readline())
	#Tgas next
	tgasheader = [f.readline() for _ in range(2)]
	tgasarr = np.zeros((npts,))
	print np.shape(tgasarr)
	for p in range(npts):
		tgasarr[p] = float(f.readline())


#####
## Plottable quantities:
# h2denarr,molabunarr,tgasarr

AU = 1.496e11

ptrad = np.sqrt(ptarr[:,0]**2+ptarr[:,1]**2)/AU
ptsph = np.sqrt(ptarr[:,0]**2+ptarr[:,1]**2+ptarr[:,2]**2)/AU
z = ptarr[:,2]/AU
absz = np.abs(ptarr[:,2])/AU

plt.figure(figsize=(24,12))


#plt.scatter(np.log10(ptrad/AU),np.log10(np.abs(ptarr[:,2])/AU),c=tgasarr,marker='.')
i = 1
if i == 0:
    plt.scatter(np.log10(ptsph),z/ptsph,c=np.log10(h2denarr),marker=',',s=4)
    plt.xlim(np.log10(6.4),np.log10(30))
    plt.ylim(-1,1)

if i ==1:
    plt.scatter(ptrad,z,c=np.log10(h2denarr),marker=',',s=1)
    plt.xlim(1,30)
    plt.ylim(-12,12)

plt.clim(10,18)
plt.colorbar()
plt.show()
