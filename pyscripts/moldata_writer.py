import numpy as np
import sys

# RUNNING THIS CODE REQUIRES TWO COMMAND LINE ARGUMENTS - v and J MAX


vmax = int(sys.argv[1]) #2
jmax = int(sys.argv[2]) #10

with open('models/co_compare/postmem/co_zz_linear_noHe.dat','r') as f:
    data_CO = [line for line in f]

# Define structure of CO file: Header 1, followed by energy levels. Header 2, followed by radiative transtions. Header 3, followed by collision table with pH2.

#Built-in assumptions: Length of hdr1 (7),hdr2 (3), hdr3 (11), hdr4(9), hdr5 (9)
#Built-in assumptions: # of coll partners (5)



hdrsum = 0
hdr1 = data_CO[:7]
hdrsum+=len(hdr1)
for l in hdr1:
    print l[:-1]
    
nlvl = int(hdr1[5])
elvls = data_CO[hdrsum:hdrsum+nlvl]

hdr2 = data_CO[hdrsum+nlvl:nlvl+hdrsum+3]
hdrsum+=len(hdr2)
for l in hdr2:
    print l[:-1]
    
ntrans = int(hdr2[1])
transns = data_CO[hdrsum+nlvl:hdrsum+nlvl+ntrans]

hdr3 = data_CO[nlvl+ntrans+hdrsum:nlvl+ntrans+hdrsum+11]
hdrsum+=len(hdr3)
for l in hdr3:
    print l[:-1]
    
nph2coll = int(hdr3[5])
ph2coll = data_CO[nlvl+ntrans+hdrsum:nlvl+ntrans+nph2coll+hdrsum]

# Continue with header 4 for oH2, header 5 for H

tmpcnt = nlvl+ntrans+nph2coll
hdr4 = data_CO[tmpcnt+hdrsum:tmpcnt+hdrsum+9]
hdrsum+=len(hdr4)
for l in hdr4:
    print l[:-1]

noh2coll = int(hdr4[3])
oh2coll = data_CO[tmpcnt+hdrsum:tmpcnt+noh2coll+hdrsum]

tmpcnt = nlvl+ntrans+nph2coll+noh2coll
hdr5 = data_CO[tmpcnt+hdrsum:tmpcnt+hdrsum+9]
hdrsum+=len(hdr5)
for l in hdr5:
    print l[:-1]

nhcoll = int(hdr5[3])
hcoll = data_CO[tmpcnt+hdrsum:tmpcnt+nhcoll+hdrsum]

newelvls = []
newtrans = []

for i in range(len(elvls)):
    tmp = elvls[i].split()
    if int(tmp[-1])<vmax+1:
        if int(tmp[-2])<jmax+1:
            newelvls.append(elvls[i])
        
for i in range(len(transns)):
    tmp = transns[i].split()
    upr = elvls[int(tmp[1])-1].split()
    lwr = elvls[int(tmp[2])-1].split()
    uv = int(upr[-1])
    uj = int(upr[-2])
    lv = int(lwr[-1])
    lj = int(lwr[-2])
    if uv<vmax+1 and lv <vmax+1 :
        if uj<jmax+1 and lj < jmax+1:
            newtrans.append(transns[i])

print len(elvls),len(newelvls)
print len(transns),len(newtrans)

newoh2coll = []
newph2coll = []
newhcoll = []

for i in range(len(ph2coll)):
    tmp = ph2coll[i].split()
    upr = elvls[int(tmp[1])-1].split()
    lwr = elvls[int(tmp[2])-1].split()
    uv = int(upr[-1])
    uj = int(upr[-2])
    lv = int(lwr[-1])
    lj = int(lwr[-2])
    if uv<vmax+1 and lv <vmax+1 :
        if uj<jmax+1 and lj < jmax+1:
            newph2coll.append(ph2coll[i])

for i in range(len(oh2coll)):
    tmp = oh2coll[i].split()
    upr = elvls[int(tmp[1])-1].split()
    lwr = elvls[int(tmp[2])-1].split()
    uv = int(upr[-1])
    uj = int(upr[-2])
    lv = int(lwr[-1])
    lj = int(lwr[-2])
    if uv<vmax+1 and lv <vmax+1 :
        if uj<jmax+1 and lj < jmax+1:
            newoh2coll.append(oh2coll[i])

for i in range(len(hcoll)):
    tmp = hcoll[i].split()
    upr = elvls[int(tmp[1])-1].split()
    lwr = elvls[int(tmp[2])-1].split()
    uv = int(upr[-1])
    uj = int(upr[-2])
    lv = int(lwr[-1])
    lj = int(lwr[-2])
    if uv<vmax+1 and lv <vmax+1 :
        if uj<jmax+1 and lj < jmax+1:
            newhcoll.append(hcoll[i])

print len(ph2coll),len(newph2coll)
print len(oh2coll),len(newoh2coll)
print len(hcoll),len(newhcoll)

# Recount all column 0 counters to account for removal of levels, transitions, collisional rates

numswp = {}
for i in range(len(newelvls)):
    numswp[int(newelvls[i][0:5])]=i+1
    newelvls[i] = str(i+1).rjust(5)+newelvls[i][5:]
    
    
    
for i in range(len(newtrans)):
    newtrans[i] = str(i+1).rjust(5)+str(numswp[int(newtrans[i][5:10])]).rjust(5)+str(numswp[int(newtrans[i][10:15])]).rjust(5)+newtrans[i][15:]
for i in range(len(newph2coll)):
    newph2coll[i] = str(i+1).rjust(6)+str(numswp[int(newph2coll[i][6:12])]).rjust(6)+str(numswp[int(newph2coll[i][12:18])]).rjust(6)+newph2coll[i][18:]
for i in range(len(newoh2coll)):
    newoh2coll[i] = str(i+1).rjust(6)+str(numswp[int(newoh2coll[i][6:12])]).rjust(6)+str(numswp[int(newoh2coll[i][12:17])]).rjust(5)+newoh2coll[i][17:]
for i in range(len(newhcoll)):
    newhcoll[i] = str(i+1).rjust(7)+str(numswp[int(newhcoll[i][7:13])]).rjust(6)+str(numswp[int(newhcoll[i][13:18])]).rjust(5)+newhcoll[i][18:]

# Recount all header numbers and overwrite

hdr1[5] = str(len(newelvls))+'\n'
hdr2[1] = str(len(newtrans))+'\n'
hdr3[5] = str(len(newph2coll))+'\n'
hdr4[3] = str(len(newoh2coll))+'\n'
hdr5[3] = str(len(newhcoll))+'\n'

writeorder = [hdr1,newelvls,hdr2,newtrans,hdr3,newph2coll,hdr4,newoh2coll,hdr5,newhcoll]
with open('co_zz_numax{0}_jmax{1}.dat'.format(vmax,jmax),'w') as f:
    for i in writeorder:
        for j in i:
            f.write(j)


