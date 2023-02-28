# Load the modules
import yt
import numpy as np
import sys

yt.enable_parallelism()

data_location = "../hdf5/ScalarFieldp_00000*"  # Data file location

# Loading dataset
ts = yt.DatasetSeries(data_location)

def _rho_dV(field, data):
    return data["chombo","cell_volume"].d*data["rho"]/data["chi"]**1.5

yt.add_field("rho_dV", _rho_dV, units="",sampling_type="local")

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):

    all_data = i.all_data()

    x_max, y_max, z_max = all_data.argmax("rho")
    loc_max = [x_max.d, y_max.d, z_max.d]
    rhosum = all_data.sum("rho_dV")

    array = [i.current_time, rhosum, loc_max]
    sto.result = array


if yt.is_root():
    timedata = []
    rhosum_data = []
    locmax_data = []

    sort = sorted(storage.items())
    for L in sort:
        timedata.append(L[1][0])
        rhosum_data.append(L[1][1])
        locmax_data.append(L[1][2])

    np.savetxt("data/time.out", timedata)
    np.savetxt("data/rhosum.out", rhosum_data)
    np.savetxt("data/locmax.out", locmax_data)

