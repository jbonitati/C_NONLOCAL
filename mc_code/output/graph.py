import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import pylab
import sys

data = pd.DataFrame.from_csv(sys.argv[1], index_col=0)
data.plot(subplots=True)
plt.ylabel('u(r)')
plt.title('Wave Functions')
plt.show()