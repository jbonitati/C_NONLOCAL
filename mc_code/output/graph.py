#Separately plots each column in a csv against the first column
#arguments: 1 - name of csv file, 2 (optional) - name of file to save plot to
#plot is automatically saved to a new file in a "figures" directory
# with the name in argument 2 or the current UTC time
# if no name is supplied

import pandas as pd
import matplotlib.pyplot as plt
import sys
import time
from shutil import copyfile

data = pd.read_csv(sys.argv[1], index_col=0)
data.plot(subplots = True,grid = True)
plt.legend(loc='best')
plt.xlabel('r (fm)')
plt.ylabel('u(r)')
plt.suptitle('Wave Functions')
if(len(sys.argv) >= 3):
  name = sys.argv[2]
  plt.savefig("figures/%s.png"%name,bbox_inches='tight')
  #copyfile("../config.ini", "inifiles/%s.ini"%name)
  print "Plot saved to figures/%s.png"%name
else:
  name = time.time()
  plt.savefig("figures/%f.png"%name, bbox_inches='tight')
  #copyfile("../config.ini", "inifiles/%f.ini"%name)
  print "Plot saved to figures/%f.png"%name
plt.show()

