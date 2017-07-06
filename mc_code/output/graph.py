#Separately plots each column in a csv against the first column
#arguments: 1 - name of csv file, 2 (optional) - name of file to save plot to
#plot is automatically saved to a new file in a "figures" directory
# with the name in argument 2 or the current UTC time
# if no name is supplied

import pandas as pd
import matplotlib.pyplot as plt
import sys
import time

data = pd.DataFrame.from_csv(sys.argv[1], index_col=0)
data.plot(subplots=True)
plt.xlabel('r (fm)')
plt.ylabel('u(r)')
plt.title('Wave Functions')
if(len(sys.argv) >= 3):
  plt.savefig("figures/%s.png"%sys.argv[2],bbox_inches='tight')
  print "Plot saved to figures/%s.png"%sys.argv[2]
else:
  plt.savefig("figures/%f.png"%time.time(), bbox_inches='tight')
  print "Plot saved to figures/%f.png"%time.time()
plt.show()

