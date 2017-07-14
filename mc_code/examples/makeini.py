#creates a template config file with specified number of channels
#arguments: 1 - name of template file to save to, 2 - number of channels

import sys

if(len(sys.argv)< 3):
  numc = input('Enter number of channels: ')
else:
  numc = sys.argv[2]
  
numc = int(numc)

#if(input('This will overwrite the config.ini in this folder. Continue? (y/n)') != "y"):
#  quit()
if(len(sys.argv) >= 2):
  name = sys.argv[1]
else:
  name = "config.ini"

f = open(name,"w+")

HEADER = ";This is a template for an input file for the multichannel R-matrix method\n\
;This template uses %d channels\n\n"%numc

SETTINGS = "[Settings] \n\
;output file will be stored in the \"output\" directory \n\
output_file = neutron_wave \n\
;note: you must have a separate section for each channel \n \
num_channels = %d \n\
entrance_channel = 1 \n\n"%numc

NUMERICAL = "[Numerical]\n\
Projectile_mass_number = 0 \n\
Target_mass_number = 0 \n\
\n\
Projectile_proton_number = 0 \n\
Target_proton_number = 0 \n\
\n\
;E (MeV) \n\
Projectile_energy = 0 \n\
\n\
Basis_size = 20\n\
Step_size = 0.01\n\
R_max = 20\n\
\n\
;a (fm)\n\
Channel_radius = 10\n\n"

LOCAL = "[local]\n\
Vv = 0\n\
rv = 0\n\
av = 1\n\
\n\
Wv = 0\n\
rwv = 0\n\
awv = 1\n\
\n\
\n\
Vd = 0\n\
rvd = 0\n\
avd = 1\n\
\n\
Wd = 0 \n\
rwd = 0\n\
awd = 1\n\
\n\
Vso = 0\n\
Rso = 0\n\
aso = 1\n\
\n\
Wso = 0\n\
Rwso = 0\n\
awso = 1\n\n"

NONLOCAL = "[Non_local]\n\
;beta = range of nonlocality\n\
beta = 0\n\
\n\
Vv = 0\n\
rv = 0\n\
av = 1\n\
\n\
Wv = 0\n\
rwv = 0\n\
awv = 1\n\
\n\
\n\
Vd = 0\n\
rvd = 0\n\
avd = 1\n\
\n\
Wd = 0\n\
rwd = 0\n\
awd = 1\n\
\n\
\n\
Vso = 0\n\
Rso = 0\n\
aso = 1\n\
\n\
Wso = 0\n\
Rwso = 0\n\
awso = 1\n\n"


f.write(HEADER)
f.write(SETTINGS)
f.write(NUMERICAL)
f.write(LOCAL)
f.write(NONLOCAL)

CHANNEL = "Angular_momentum = 0\n\
Total_angular_momentum = 0\n\
Energy = 0\n\n"

for i in range(numc):
  f.write("[Channel%d]\n"%(i+1))
  f.write(CHANNEL)
  
COUPLING = "beta = 0\n\
V = 0\n\
r = 0\n\
a = 1\n\n"
  
for i in range(numc):
  j = i+2
  while(j <= numc):
    f.write("[coupling%d%d]\n"%((i+1),(j)))
    f.write(COUPLING)
    j = j+1
