import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

# ============================================================
# READ FORTRAN FILE
# ============================================================

f = FortranFile('tl.grid', 'r')

records = []

while True:

    try:

        rec = f.read_reals(np.float32)

        records.append(rec)

    except:
        break

print('records =', len(records))

# ============================================================
# REMOVE HEADER
# ============================================================

records = records[1:]

# ============================================================
# STACK INTO ARRAY
# ============================================================

TL = np.vstack(records)

print('TL shape =', TL.shape)

# ============================================================
# TRANSPOSE
# ============================================================

TL = TL.T

print('Transposed shape =', TL.shape)

# ============================================================
# BUILD AXES
# ============================================================

nr = TL.shape[1]
nz = TL.shape[0]

r = np.linspace(0,50,nr)      # km
z = np.linspace(0,4000,nz)   # m

# ============================================================
# PLOT
# ============================================================

plt.figure(figsize=(12,6))

plt.imshow(
    TL,
    extent=[r[0],r[-1],z[-1],z[0]],
    aspect='auto',
    cmap='jet'
)

plt.xlabel('Range (km)')
plt.ylabel('Depth (m)')

plt.title('RAM Transmission Loss')

plt.colorbar(label='TL')

plt.show()
