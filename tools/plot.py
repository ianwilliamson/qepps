import argparse
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

parser=argparse.ArgumentParser(description='QEPPS Plot Utility')
parser.add_argument('files',nargs='+',help='Data text files to plot')
parser.add_argument('-out',default='',help='Save filename')
parser.add_argument('-col',type=int,default=1,help='number of data columns to plot')
args=parser.parse_args()

ax = plt.gca()
for fn in args.files:
  data=np.loadtxt(fn,comments='#',delimiter=', ',dtype=np.complexfloating)
  freqs=np.real(data[:,0])
  neff =data[:,1:args.col+1]
  plt.plot( freqs,  np.real(neff), linewidth=2, )
  plt.plot( freqs, -np.imag(neff), linewidth=2, linestyle='--', color=ax.lines[-1].get_color() )

plt.xlabel('Frequency')
plt.ylabel('n')

if args.out:
  plt.savefig(args.out)
else:
  plt.show()

