import numpy as np
from matplotlib import pyplot as plt
import scipy.interpolate
from matplotlib.colors import Normalize
from matplotlib import cm

### FUNCTION DEFINITIONS ###############################################

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, theta):
    x = rho * np.sin(theta)
    z = rho * np.cos(theta)
    return(x, z)

def sph2cart(rho, theta, phi):
    x = rho * np.cos(phi) * np.sin(theta)
    y = rho * np.sin(phi) * np.sin(theta)
    z = rho * np.cos(theta)
    return(np.array([x, y, z]))

def sph2cartvec(vec,tht,phi):
    #Vec given in [r,tht,phi] order. want [phi,tht,r]
    #Output in x,y,z order
    vec = vec[::-1]
    transform = np.array([
        [-np.sin(phi), -np.sin(tht)*np.cos(phi), np.cos(tht)*np.cos(phi)],
        [np.cos(phi), -np.sin(tht)*np.sin(phi), np.cos(tht)*np.sin(phi)],
        [0., np.cos(tht), np.sin(tht)]])
    return np.inner(transform,vec)

def linefind(r, t, tmax):
    return (t+r*tmax)


#### LOAD RADLITE GRID INFO #############################################

radlite_home = '/Users/tpauly/projects_unsync/radall/'
model = raw_input('Which radlite model directory?')
radlite_model_dir = radlite_home+model+'/'

with open(radlite_model_dir+'radius.inp','r') as f:
    rad = [line for line in f]
    rad = np.array([float(x) for x in rad[2:]])
   
with open(radlite_model_dir+'theta.inp','r') as g:
    tht = [line for line in g]
    tht = np.array([float(y) for y in tht[2:]])

nrad = len(rad)
ntht = len(tht)

with open(radlite_model_dir+'density.inp','r') as h:
    den = [line for line in h]
    den = np.array([float(y) for y in den[1:]])

with open(radlite_model_dir+'abundance.inp','r') as j:
    abn = [line for line in j]
    abn = np.array([float(y.split()[0]) for y in abn[1:]])
    
with open(radlite_model_dir+'temperature.inp','r') as k:
    tem = [line for line in k]
    tem = np.array([float(y) for y in tem[1:]])
    
with open(radlite_model_dir+'dusttemp_final.dat','r') as l:
    dtem = [line for line in l]
    dtem = np.array([float(y) for y in dtem[3:]])
    
with open(radlite_model_dir+'turbulence.inp','r') as m:
    trb = [line for line in m]
    trb = np.array([float(y) for y in trb[2:]])
    
with open(radlite_model_dir+'velocity.inp','r') as n:
    temp = [line for line in n]
    vel = np.zeros((len(temp)-1,3))
    for i in range(1,len(temp)):
        vel[i-1] = np.array([float(x) for x in temp[i].split()])
        
with open('najita_hydrogen_r0p25au.txt','r') as f:
    hydinnr = np.loadtxt(f)

with open('najita_hydrogen_r1au.txt','r') as f:
    hydmidr = np.loadtxt(f)

with open('najita_hydrogen_r20au.txt','r') as f:
    hydoutr = np.loadtxt(f)


#### Need to construct a matrix of vertical column density to use Najita H/H2 figures. 
#### Assume that majority of column is found in small angles, such that taking 
#### r=const, vary theta gives ~vertical column.

den2d = den.reshape((nrad,ntht))/1.674e-24
hden2d = np.zeros_like(den2d)
h2den2d = np.zeros_like(den2d)

r = np.zeros_like(den2d)
z = np.zeros_like(den2d)
cum_vertcol = np.zeros_like(den2d)
for i in range(nrad):
    for j in range(ntht):
        r[i,j] = rad[i]*np.sin(tht[j])
        z[i,j] = rad[i]*np.cos(tht[j])
        if j == 0:
            #First segment of vertical path starts at tht != 0,
            #so first segment goes from r=0, z=rad to r_ij, z_ij.
            #Use density of d_ij for this path
            seg_col = den2d[i,j]*np.sqrt(r[i,j]**2+(rad[i]-z[i,j])**2)
            cum_vertcol[i,j] = seg_col
        else:
            seg_col = den2d[i,j]*np.sqrt((r[i,j]-r[i,j-1])**2+(z[i,j]-z[i,j-1])**2)
            cum_vertcol[i,j] = cum_vertcol[i,j-1]+seg_col

def linefit(arr, loc, step, idx):
    if step > 0:
        slope = (arr[loc+step,idx]-arr[loc,idx])/(arr[loc+step,0]-arr[loc,0])
        intercept = arr[loc,idx]-(slope*arr[loc,0])
    else:
        slope = (arr[loc,idx]-arr[loc+step,idx])/(arr[loc,0]-arr[loc+step,0])
        intercept = arr[loc,idx]-(slope*arr[loc,0])
    return slope,intercept

def sanity(arg,uplim):
    if arg > uplim:
        return uplim
    elif arg < 1.e-10:
        return 1.e-10
    else:
        return arg

for i in range(len(rad)):
    for j in range(len(tht)):
        col = cum_vertcol[i,j]
        if rad[i] < 1.496e13/4.:
            #Just use R=0.25 AU Najita line
            mincol = np.argmin((hydinnr[:,0]-col)**2)
            st = -1 if (hydinnr[mincol,0]-col > 0) else 1
            if mincol == len(hydinnr) - 1:
                st = -1
            sl1,b1 = linefit(hydinnr,mincol,st,1)
            sl2,b2 = linefit(hydinnr,mincol,st,2)
            hden2d[i,j] = sanity(b1 + sl1*col,1.)
            h2den2d[i,j] = sanity(b2 + sl2*col,0.5)
        elif rad[i] < 1.496e13:
            #Interp between R=0.25 and R=1 AU Najita data
            mincol1 = np.argmin((hydinnr[:,0]-col)**2)
            st1 = -1 if (hydinnr[mincol1,0]-col > 0) else 1
            if mincol1 == len(hydinnr) - 1:
                st1 = -1
            insl1,inb1 = linefit(hydinnr,mincol1,st1,1)
            insl2,inb2 = linefit(hydinnr,mincol1,st1,2)
            inh = inb1 + insl1*col
            inh2 = inb2 + insl2*col
            ####
            mincol2 = np.argmin((hydmidr[:,0]-col)**2)
            st2 = -1 if (hydmidr[mincol2,0]-col > 0) else 1
            if mincol2 == len(hydmidr) - 1:
                st2 = -1

            mdsl1,mdb1 = linefit(hydmidr,mincol2,st2,1)
            mdsl2,mdb2 = linefit(hydmidr,mincol2,st2,2)
            mdh = mdb1 + mdsl1*col
            mdh2 = mdb2 + mdsl2*col
            ####
            hsl = (mdh-inh)/(1.496e13-(1.496e13/4.))
            hb = inh-hsl*1.496e13/4.
            h2sl = (mdh2-inh2)/(1.496e13-(1.496e13/4.))
            h2b = inh2-h2sl*1.496e13/4.
            hden2d[i,j] = sanity(hb + hsl*rad[i],1.)
            h2den2d[i,j] = sanity(h2b + h2sl*rad[i],0.5)
        elif rad[i] < 20.*1.496e13:
            #Interp between R=1 and R=20 AU Najita data
            mincol1 = np.argmin((hydmidr[:,0]-col)**2)
            st1 = -1 if (hydmidr[mincol1,0]-col > 0) else 1
            if mincol1 == len(hydmidr) - 1:
                st1 = -1

            mdsl1,mdb1 = linefit(hydmidr,mincol1,st1,1)
            mdsl2,mdb2 = linefit(hydmidr,mincol1,st1,2)
            mdh = mdb1 + mdsl1*col
            mdh2 = mdb2 + mdsl2*col
            ####
            mincol2 = np.argmin((hydoutr[:,0]-col)**2)
            st2 = -1 if (hydoutr[mincol2,0]-col > 0) else 1
            if mincol2 == len(hydoutr) - 1:
                st2 = -1

            ousl1,oub1 = linefit(hydoutr,mincol2,st2,1)
            ousl2,oub2 = linefit(hydoutr,mincol2,st2,2)
            ouh = oub1 + ousl1*col
            ouh2 = oub2 + ousl2*col
            ####
            hsl = (ouh-mdh)/(20.*1.496e13-(1.496e13))
            hb = mdh-hsl*1.496e13
            h2sl = (ouh2-mdh2)/(20.*1.496e13-(1.496e13))
            h2b = mdh2-h2sl*1.496e13
            hden2d[i,j] = sanity(hb + hsl*rad[i],1.)
            h2den2d[i,j] = sanity(h2b + h2sl*rad[i],0.5)
        else:
            #Just use R=20 AU Najita line
            mincol = np.argmin((hydoutr[:,0]-col)**2)
            st = -1 if (hydoutr[mincol,0]-col > 0) else 1
            if mincol == len(hydoutr) - 1:
                st = -1

            sl1,b1 = linefit(hydoutr,mincol,st,1)
            sl2,b2 = linefit(hydoutr,mincol,st,2)
            hden2d[i,j] = sanity(b1 + sl1*col,1.)
            h2den2d[i,j] = sanity(b2 + sl2*col,0.5)

###### NAME OF OUTPUT FILE  ################################################
outfile = open('radlite_grid'+model+'.dat','w')


print np.min(h2den2d),np.max(h2den2d)

phi = [0.0]
nphi = len(phi)
for i in range(len(rad)):
    #Add loop to generate null inner point for gap creation
    if i==0:
        for j in range(len(tht)):
            tmprad = rad[0]*(rad[0]/rad[1])
            x, y, z = sph2cart(tmprad, tht[j], phi[0])/100.  #Div by 100 to go from cm to m
            outstr1 = "{:+.5E}  {:+.5E}  {:+.5E}  ".format(x, y, z)
            ln = linefind(i,j,ntht)
            outstr2 = "{:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}  ".format(1.e-24,abn[ln],tem[ln],dtem[ln],trb[ln]*1.e3)
            velx,vely,velz = sph2cartvec(vel[ln],tht[j],phi[0])/100.  #Div by 100 to go from cm/s to m/s
            #Added habun,h2abun to outstr3
            outstr3 = "{:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}\n".format(velx,vely,velz,hden2d[i,j],h2den2d[i,j])
            outfile.write(outstr1+outstr2+outstr3)
            x, y, z = sph2cart(tmprad, (np.pi) - tht[j], phi[0])/100.
            outstr4 = "{:+.5E}  {:+.5E}  {:+.5E}  ".format(x, y, z)
            velx,vely,velz = sph2cartvec(vel[ln],(np.pi) - tht[j],phi[0])/100.
            #Added habun,h2abun to outstr5
            outstr5 = "{:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}\n".format(velx,vely,velz,hden2d[i,j],h2den2d[i,j])
            #outfile.write(outstr4+outstr2+outstr5)


    for j in range(len(tht)):
        x, y, z = sph2cart(rad[i], tht[j], phi[0])/100.  #Div by 100 to go from cm to m
        outstr1 = "{:+.5E}  {:+.5E}  {:+.5E}  ".format(x, y, z)
        ln = linefind(i,j,ntht)
        outstr2 = "{:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}  ".format(den[ln],abn[ln],tem[ln],dtem[ln],trb[ln]*1.e3)
        velx,vely,velz = sph2cartvec(vel[ln],tht[j],phi[0])/100.  #Div by 100 to go from cm/s to m/s
        outstr3 = "{:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}\n".format(velx,vely,velz,hden2d[i,j],h2den2d[i,j])
        outfile.write(outstr1+outstr2+outstr3)
        x, y, z = sph2cart(rad[i], (np.pi) - tht[j], phi[0])/100.
        outstr4 = "{:+.5E}  {:+.5E}  {:+.5E}  ".format(x, y, z)
        velx,vely,velz = sph2cartvec(vel[ln],(np.pi) - tht[j],phi[0])/100.
        outstr5 = "{:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}  {:+.5E}\n".format(velx,vely,velz,hden2d[i,j],h2den2d[i,j])
        #outfile.write(outstr4+outstr2+outstr5)

outfile.close()

fullrad = np.insert(rad,0,tmprad,axis=0)
fullrad = fullrad/100.
with open('radius'+model+'.dat','w') as f:
	for i in range(len(fullrad)):
		f.write('{:e}\n'.format(fullrad[i]))

with open('theta'+model+'.dat','w') as f:
	for i in range(len(tht)):
		f.write('{:e}\n'.format(tht[i]))