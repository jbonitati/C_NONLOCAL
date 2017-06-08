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


data=loadtxt("example1.txt",float)

data2=loadtxt("example2.txt",float)

#print(data)
#data=loadtxt("misc3.txt",float)
scale_pow = 2
x0=data[:,0]
y0=data[:,1]
z0=data[:,2]
x1=data2[:,0]
y1=data2[:,1]
z1=data2[:,2]
plt.suptitle('Ca40 scattering state binding energy=5 MeV l=0', fontsize=9, fontweight='bold')


plt.plot(x0,y0,'*',label='old real')
plt.legend(loc='lower left')
plt.plot(x1,y1,'o',label='new real')
plt.legend(loc='lower left')
plt.plot(x0,z0,'*',label='old imag')
plt.legend(loc='lower left')
plt.plot(x1,z1,'o',label='new imag')
plt.legend(loc='lower left')
plt.ylabel('log(\phi(r))')
plt.xlabel('r[fm]')
plt.show()


