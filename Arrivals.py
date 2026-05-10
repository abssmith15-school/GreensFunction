# ============================================================
# Kraken Broadband Time Front Synthesis
# ============================================================

from os import system
import numpy as np
import matplotlib.pyplot as plt

from wkrakenenvfil import *
from read_shd import *

# ============================================================
# BASIC PARAMETERS
# ============================================================

fc = 100.0
Q  = 10.0
bw = fc / Q

rmax = 20e3                 # 20 km
Dmax = 150.0                # shallow water

zs = np.array([50.0])       # SOURCE DEPTH = 50 m

# ============================================================
# TIME / FFT SETUP
# ============================================================

c0 = 1500.0

T  = 5.0
fs = 4 * fc

dt = 1/fs
N  = int(fs * T)

df = 0.2

frq = np.arange(df, bw, df)
frq = np.concatenate((-frq[::-1], [0], frq)) + fc

nf = len(frq)

wind = np.sinc((frq - fc)/bw)

# ============================================================
# RECEIVER GRID
# ============================================================

zarray = np.arange(1, Dmax+1, 1)     # RECEIVER DEPTHS
nza = len(zarray)

rarray = np.linspace(1, rmax/1000, 200)  # km

# ============================================================
# STORAGE
# ============================================================

psif = np.zeros((nza, nf), dtype=complex)

# ============================================================
# SURFACE DATA
# ============================================================

bc = 'V'
properties = []
reflection = []

surface_data = {
    "bc":bc,
    "properties":properties,
    "reflection":reflection
}

# ============================================================
# SCATTER DATA
# ============================================================

bumden = []
eta = []
xi = []

scatter_data = {
    "bumden":bumden,
    "eta":eta,
    "xi":xi
}

# ============================================================
# SSP DATA
# ============================================================

sspdata = np.loadtxt("elba.ssp")

z = sspdata[:,0]
c = sspdata[:,1]

nz = z.size

cs = np.zeros(nz)
rho = np.ones(nz)
apt = cs
ast = cs

type  = 'H'
itype = 'N'

nmesh = 0
sigma = 0.0

clow  = 1450.0
chigh = 2500.0

cdata = np.array([z,c,cs,rho,apt,ast])

ssp_data = {
    "cdata":cdata,
    "type":type,
    "itype":itype,
    "nmesh":nmesh,
    "sigma":sigma,
    "clow":clow,
    "chigh":chigh,
    "zbottom":max(z)
}

# ============================================================
# BOTTOM DATA
# ============================================================

layer_info = np.array([
    [98.0,1600.0,130.0,1.49,0.1,1.70],
    [103.0,1800.0,500.0,1.90,0.90,2.50],
    [128.0,2500.0,900.0,2.40,0.01,0.01]
])

layerp = np.array([
    [nmesh,0.0,layer_info[1,0]],
    [nmesh,0.0,layer_info[2,0]]
])

layert = 'HH'

properties = np.array(layer_info[2,:])

m1 = np.array([layer_info[0,:],layer_info[0,:]])
m2 = np.array([layer_info[1,:],layer_info[1,:]])

bdata = np.array([m1,m2])

bdata[0,1,0] = layer_info[1,0]
bdata[1,1,0] = layer_info[2,0]

units = 'W'
bc = 'A'
sigma = 0.0

bottom_data = {
    "n":3,
    "layerp":layerp,
    "layert":layert,
    "properties":properties,
    "bdata":bdata,
    "units":units,
    "bc":bc,
    "sigma":sigma
}

# ============================================================
# FIELD DATA
# ============================================================

rp = 0
np_ = 1

m = 999

rmodes = 'A'
stype  = 'R'
thorpe = 'T'
finder = ' '

dr = np.zeros(nza)

field_data = {
    "rmax":rmax/1000,
    "nrr":len(rarray),
    "rr":rarray,
    "rp":rp,
    "np":np_,
    "m":m,
    "rmodes":rmodes,
    "stype":stype,
    "thorpe":thorpe,
    "finder":finder,
    "rd":zarray,
    "dr":dr,
    "nrd":nza
}

# ============================================================
# FREQUENCY LOOP
# ============================================================

for i, freq in enumerate(frq):

    print(f'\nRunning frequency {freq:.2f} Hz')

    # --------------------------------------------------------
    # SOURCE DATA
    # --------------------------------------------------------

    source_data = {
        "zs": zs,
        "f": freq
    }

    # --------------------------------------------------------
    # WRITE ENV FILE
    # --------------------------------------------------------

    wkrakenenvfil(
        'elba',
        'Broadband Elba',
        source_data,
        surface_data,
        scatter_data,
        ssp_data,
        bottom_data,
        field_data
    )

    # --------------------------------------------------------
    # WRITE FLP FILE
    # --------------------------------------------------------

    with open('elba.flp', 'w', newline='\n') as flp:

        # Title
        flp.write("'elba'\n")

        # Run type
        flp.write("'RA'\n")

        # Frequency
        flp.write(f"{freq:.0f}\n")

        # NPROF
        flp.write("1\n")

        # Profile range
        flp.write("0\n")

        # Number of receiver ranges
        flp.write(f"{len(rarray)}\n")

        # Receiver ranges (km)
        flp.write(f"1 {rmax/1000:.0f} /\n")

        # Number of source depths
        flp.write(f"{len(zs)}\n")

        # Source depth
        flp.write(f"{zs[0]:.0f}\n")

        # Number of receiver depths
        flp.write(f"{len(zarray)}\n")

        # Receiver depth limits
        flp.write(f"1 {Dmax:.0f} /\n")

        # Number of receiver range offsets
        flp.write(f"{len(zarray)}\n")

        # Receiver range offsets
        flp.write("0 /\n")

    # --------------------------------------------------------
    # RUN KRAKEN
    # --------------------------------------------------------

    ret1 = system("./krakenc.exe elba")
    print("kraken return =", ret1)

    # --------------------------------------------------------
    # RUN FIELD
    # --------------------------------------------------------

    ret2 = system("./field.exe elba")
    print("field return =", ret2)

    # --------------------------------------------------------
    # READ SHD FILE
    # --------------------------------------------------------

    PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure = read_shd('elba')

    print("pressure shape =", np.shape(pressure))

    p = pressure[0,0,:,:]

    print("max pressure =", np.nanmax(np.abs(p)))

    # --------------------------------------------------------
    # LAST RANGE
    # --------------------------------------------------------

    pfield = p[:, -1]

    psif[:, i] = pfield

# ============================================================
# SYNTHESIZE ARRIVALS
# ============================================================

tdelay = rmax / c0

ptz = np.zeros((nza, N), dtype=complex)

nyqst = int(np.ceil((nf + 1)/2))

for iz in range(nza):

    data = (
        wind
        * np.conj(psif[iz,:])
        * np.exp(1j * 2*np.pi * frq * tdelay)
    )

    data = np.concatenate((
        data[nyqst:nf],
        np.zeros(N - nf),
        data[0:nyqst]
    ))

    ptz[iz,:] = np.fft.ifft(data)

# ============================================================
# TIME AXIS
# ============================================================

taxis = tdelay + np.arange(N)/fs

# ============================================================
# PLOT ARRIVALS
# ============================================================

pk = np.max(np.abs(ptz))

if pk <= 0:
    raise ValueError("Synthesized field is zero.")

data_db = 20*np.log10(
    np.maximum(np.abs(ptz)/pk, 1e-30)
)

data_db[data_db < -70] = -70

extent = [
    taxis[0],
    taxis[-1],
    zarray[-1],
    zarray[0]
]

plt.figure(figsize=(10,8))

plt.imshow(
    data_db,
    extent=extent,
    aspect='auto',
    cmap='jet',
    vmin=-70,
    vmax=0
)

plt.colorbar(label='Relative Level (dB)')

plt.xlabel('Travel Time (s)')
plt.ylabel('Depth (m)')
plt.title('Kraken Broadband Time Front')

plt.show()
