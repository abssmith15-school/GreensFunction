#==================================================================
#
#  KRAKEN: Elba waveguide (single bottom layer)
#
#==================================================================

from os import *
import sys
from numpy import *
from matplotlib.pyplot import *

from wkrakenenvfil import *
from read_shd import *
from read_modes import *

print('Elba waveguide:')

case_title = 'Elba waveguide'

#==================================================================
#
# Basic parameters
#
#==================================================================

freq = 100.0       # Frequency [Hz]

Dmax = 350.0       # Water depth [m]

cw   = 1500.0      # Water sound speed
cb   = 1700.0      # Bottom compressional speed

rhow = 1.0         # Water density
rhob = 1.8         # Bottom density

w = 2*pi*freq

#==================================================================
#
# Source
#
#==================================================================

zs = array([50.0])   # Source depth
rs = array([0.0])

source_data = {
    "zs": zs,
    "f": freq
}

#==================================================================
#
# Computational grid
#
#==================================================================

rmax = 20000.0
rmaxkm = rmax / 1000

nra = 201

rarray = linspace(0, rmax, nra)
rarraykm = rarray / 1000

zarray = arange(0, Dmax + 1, 1)
nza = zarray.size

#==================================================================
#
# Surface data
#
#==================================================================

bc = 'V'

properties = []
reflection = []

surface_data = {
    "bc": bc,
    "properties": properties,
    "reflection": reflection
}

#==================================================================
#
# Scatter data
#
#==================================================================

bumden = []
eta    = []
xi     = []

scatter_data = {
    "bumden": bumden,
    "eta": eta,
    "xi": xi
}

#==================================================================
#
# Sound speed profile
#
#==================================================================

sspdata = loadtxt("elba.ssp")

z = sspdata[:,0]
c = sspdata[:,1]

zmax = max(z)

# Save original SSP
z0 = z
c0 = c

nz = z.size

# Water is acoustic only
cs  = zeros(nz)

rho = ones(nz)

apt = zeros(nz)
ast = zeros(nz)

type  = 'H'
itype = 'N'

# Mesh points
nmesh = 0

sigma = 0.0

clow  = 1450.0
chigh = 1800.0

cdata = array([z, c, cs, rho, apt, ast])

ssp_data = {
    "cdata": cdata,
    "type": type,
    "itype": itype,
    "nmesh": nmesh,
    "sigma": sigma,
    "clow": clow,
    "chigh": chigh,
    "zbottom": zmax
}

#==================================================================
#
# SINGLE BOTTOM LAYER
#
# [ z_top , cp , cs , rho , ap , as ]
#
#==================================================================

layer_info = array([
    [300.0, 1650.0, 0.0, 1.8, 0.3, 0.0]
])

#==================================================================
#
# Bottom data
#
#==================================================================

layerp = array([
    [nmesh, 0.0, layer_info[0,0]]
])

layert = 'H'

properties = array(layer_info[0,:])

m1 = array([
    layer_info[0,:],
    layer_info[0,:]
])

bdata = array([m1])

units = 'W'

bc = 'A'

sigma = 0.0

bottom_data = {
    "n": 1,
    "layerp": layerp,
    "layert": layert,
    "properties": properties,
    "bdata": bdata,
    "units": units,
    "bc": bc,
    "sigma": sigma
}

#==================================================================
#
# Field data
#
#==================================================================

rp     = 0
np     = 1

m      = 999

rmodes = 'A'

stype  = 'R'

thorpe = 'T'

finder = ' '

dr = zeros(nza)

field_data = {
    "rmax": rmaxkm,
    "nrr": nra,
    "rr": rarraykm,
    "rp": rp,
    "np": np,
    "m": m,
    "rmodes": rmodes,
    "stype": stype,
    "thorpe": thorpe,
    "finder": finder,
    "rd": zarray,
    "dr": dr,
    "nrd": nza
}

#==================================================================
#
# Write environment
#
#==================================================================

print("Writing environmental file...")

wkrakenenvfil(
    'elba',
    case_title,
    source_data,
    surface_data,
    scatter_data,
    ssp_data,
    bottom_data,
    field_data
)

#==================================================================
#
# Run KRAKEN
#
#==================================================================

print("Running KRAKEN...")

system("./krakenc.exe elba")

system("cp field.flp elba.flp")

system("./field.exe elba < elba.flp")

#==================================================================
#
# Read output
#
#==================================================================

print("Reading output data...")

PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure = read_shd('elba')

Modes = read_modes('elba.mod', freq)

k   = Modes['k']
phi = Modes['phi']
z   = Modes['z']

#==================================================================
#
# Transmission loss
#
#==================================================================

p = squeeze(pressure, axis=(0,1))

tl = -20 * log10(abs(p))

#==================================================================
#
# Plot TL
#
#==================================================================

figure(1)

imshow(
    tl,
    extent=[0, rmax, Dmax, 0],
    aspect='auto',
    cmap='jet_r',
    vmin=40,
    vmax=80
)

colorbar()



xlabel('Range (m)')

ylabel('Depth (m)')

title('KRAKEN - Elba waveguide')

ylim(Dmax, 0)

#==================================================================
#
# Plot first 4 modes
#
#==================================================================

figure(2)

nmodes_plot = min(4, phi.shape[1])

for i in range(nmodes_plot):

    rphi = real(phi[:,i])

    iphi = imag(phi[:,i])

    thetitle = 'Mode ' + str(i+1)

    subplot(1, nmodes_plot, i+1)

    plot(rphi, z, 'k')
    plot(iphi, z, 'r--')

    title(thetitle)

    ylim(Dmax, 0)

    grid(True)

subplot(1, nmodes_plot, 1)

ylabel('Depth (m)')

show()

print("done.")
