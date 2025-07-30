import h5py
import numpy as np

### Writing to HDF5 files ###

def _save_dict_to_h5(group, data_dict):
    for key, value in data_dict.items():
        if isinstance(value, dict):
            # If it's a dictionary, create a new group
            new_group = group.create_group(key)
            _save_dict_to_h5(new_group, value) # Recursively call for nested dict
        else:
            # # Check if numpy array or standard python type
            # if not isinstance(value, np.ndarray) or not isinstance(value, (int, float, str, bool)):
            #     raise ValueError(f"Value for key '{key}' is not a numpy array.")
            # if value.ndim > 2:
            #     raise ValueError(f"Value for key '{key}' has more than 2 dimensions, which is not supported.")
            # If it's an array, create a dataset
            group.create_dataset(key, data=value)


def save_to_hdf5(my_data, filename):
    """
    Save a dictionary of numpy arrays to an HDF5 file.

    Parameters:
    my_data (dict): Dictionary where keys are dataset names and values are numpy arrays.
    file_name (str): Name of the output HDF5 file.
    """
    # 1. Open a new HDF5 file
    # 2. Create datasets for each array in the dictionary
    with h5py.File(filename, 'w') as f:
        _save_dict_to_h5(f, my_data)


### Reading from HDF5 files ###

def _load_dict_from_h5(group):
    data_dict = {}
    for key in group.keys():
        item = group[key]
        if isinstance(item, h5py.Group):
            # If it's a group, recursively load it
            data_dict[key] = _load_dict_from_h5(item)
        else:
            # Otherwise, it's a dataset
            data = item[()]
            if isinstance(data, bytes):
                # Decode bytes to string if necessary
                data_dict[key] = data.decode("utf-8")
            else:
                data_dict[key] = data
    return data_dict


def load_from_hdf5(filename):
    """
    Recursively load a nested dictionary of numpy arrays from an HDF5 file.

    Parameters:
    filename (str): Name of the HDF5 file to read.
    Returns:
    dict: A dictionary where keys are dataset names and values are numpy arrays.
    """

    with h5py.File(filename, 'r') as f:

        return _load_dict_from_h5(f)
    

def load_key_from_hdf5(filename, key):
    """
    Load a specific key from an HDF5 file.

    Parameters:
    filename (str): Name of the HDF5 file to read.
    key (str): The key to load from the file.

    Returns:
    The data associated with the specified key.
    """
    with h5py.File(filename, 'r') as f:
        if key in f:
            if isinstance(f[key], h5py.Group):
                # If it's a group, recursively load it
                return _load_dict_from_h5(f[key])
            else:
                return f[key][()]
        else:
            raise KeyError(f"Key '{key}' not found in the file '{filename}'.")

