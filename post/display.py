import numpy as np
import matplotlib.pyplot as plt
import argparse as ap

parser = ap.ArgumentParser( description="Display a component of the composition vs. time." )
parser.add_argument( "-fn", "--filename", 
                     help="Filename of output data to read (default: pasr.out).", 
                     default="pasr.out" ) 
parser.add_argument( "-var", "--variable",
                     help="Variable to plot as it appears in header of output file (default: T).",
                     default="T")
args = parser.parse_args()

# read header line
hdrstr = open( args.filename ).readline()
header = hdrstr.split()
fcn_ind = header.index( args.variable )

# read points
profile = np.loadtxt( args.filename, skiprows=1 )
t = profile[:,0]
fcn = profile[:,fcn_ind]

plt.plot( t, fcn, lw=3 )

# label figure
plt.title( args.variable + ' vs. Time', fontsize=18, fontweight='bold')
plt.xlabel( 't (s)', fontsize=16, fontweight='bold' )
plt.ylabel( args.variable, fontsize=16, fontweight='bold' )
plt.legend( loc='upper left' )

# plot
plt.box('on')
plt.grid('on')
plt.show()
