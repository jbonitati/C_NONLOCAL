import numpy as np
import numpy
import scipy
import math
import cmath
from scipy import special
from scipy.special import legendre
import scipy.constants as sc
from scipy.integrate import quad
from scipy.integrate import quadrature
from scipy.misc import derivative

from scipy.special import spherical_jn 
from scipy import integrate
from scipy.special import sph_jn 
from scipy.special import spherical_yn
from numpy.linalg import inv
from numpy.linalg import pinv
from numpy.linalg import tensorinv
from matplotlib.legend_handler import HandlerLine2D
import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt
#from le import GaussLegendreWeights
import cmath



from numpy import loadtxt
import matplotlib.pyplot as plt
from matplotlib import pyplot
import matplotlib.ticker as mtick
import numpy as np
from matplotlib.legend_handler import HandlerLine2D


data=loadtxt("neutron_distorted_wave.txt",float)

#print(data)
#data=loadtxt("misc3.txt",float)
scale_pow = 2
x0=data[:,0]
y01=data[:,1]
y011=[]
for i in y01:
   y011.append(i**2.0)

x1=data[:,0]
y11=data[:,2]
y022=[]
for i in y11:
    y022.append(i**2.0)
from operator import add

y0= map(add, y011, y022)


plt.suptitle('Ca40 scattering state binding energy=5 MeV l=0', fontsize=9, fontweight='bold')


plt.plot(x0,y0,'*',label='weichuan')
plt.legend(loc='lower left')

  
plt.ylabel('log(\phi(r))')
plt.xlabel('r[fm]')
plt.show()


