import h5py
import numpy as np
import yt
from typing import Tuple
import matplotlib.pyplot as plt

# Open the HDF5 file in read-only mode

base_file_name = "/work/ta084/ta084/harveymorris/single_monopole_flat_test/hdf5/ScalarFieldp_000" #000.3d.hdf5"

def get_centre_of_monopole(rho: np.ndarray) -> Tuple[float, float]:
    total_sum = np.sum(rho)
    idx1_mean: float = 0
    idx2_mean: float = 0
    for idx1 in range(rho.shape[0]):
        for idx2 in range(rho.shape[1]):
            idx1_mean += rho[idx1, idx2] * idx1 / total_sum
            idx2_mean += rho[idx1, idx2] * idx2 / total_sum
    return idx1_mean, idx2_mean

def pad_string_number(num: int) -> str:
    number = str(num)
    while len(number) < 3:
        number = '0' + number
    return number


def get_z_location(file_name: str) -> float:
    ds = yt.load(file_name)
    monopole_rho = ds.slice(0, 256)["rho"].to_ndarray().reshape(128, 128)#[:, :64]
    np.save('monopole_rho1.npy', monopole_rho)
    print(get_centre_of_monopole(rho=monopole_rho)[1])
    

get_z_location(file_name=base_file_name + "000.3d.hdf5")




# Convert the dataset to a numpy array
#data = dataset.to_ndarray()

# Print the shape of the array
#print(data.shape)