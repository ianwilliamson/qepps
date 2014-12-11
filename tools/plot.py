import argparse
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

parser=argparse.ArgumentParser(description='QEPPS Plot Utility')
parser.add_argument('files',nargs='+',help='QEEPS output text file(s)')
parser.add_argument('-o','--output',default='',help='Save filename, if not specified will be shown interactively')
parser.add_argument('--imag',action='store_true',help='Flag specifying if the imaginary component of the eigenvalues should be plotted')
parser.add_argument('-c','--columns',type=int,default=1,help='Number of data columns to be plotted')
args=parser.parse_args()

plt.style.use('https://enw.me/pyplot/ian.mplstyle')
ax = plt.gca()
for fn in args.files:
  data=np.loadtxt(fn,comments='#',delimiter=',',dtype=complex)
  freqs=np.real(data[:,0])/1E12
  neff =data[:,1:args.columns+1]
  plt.plot( freqs,  np.abs(np.real(neff)))
  if args.imag:
    plt.plot( freqs, -np.imag(neff), linestyle='--',color=ax.lines[-1].get_color() )
  
  fmin=np.amin(freqs)
  fmax=np.amax(freqs)
  plt.xlim( (fmin,fmax) ) 
  

plt.xlabel('Frequency (THz)')
plt.ylabel('k/k_0')

plt.legend

if args.output:
  plt.savefig(args.output)
else:
  plt.show()

