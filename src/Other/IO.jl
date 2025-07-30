using HDF5

### Writing to HDF5 files ###

function _save_dict_to_h5(group, data_dict)
    for (key, value) in data_dict
        if isa(value, Dict)
            # If it's a dictionary, create a new group
            new_group = create_group(group, string(key))
            _save_dict_to_h5(new_group, value)  # Recursively call for nested dict
        else
            # If it's an array or scalar, create a dataset
            if isa(value, Symbol)
                # Convert Symbol to String for HDF5 compatibility
                group[string(key)] = string(value)
            else
                group[string(key)] = value
            end
        end
    end
end

function save_to_hdf5(my_data, filename)
    """
    Save a dictionary of arrays to an HDF5 file.

    Parameters:
    my_data (Dict): Dictionary where keys are dataset names and values are arrays.
    filename (String): Name of the output HDF5 file.
    """
    h5open(filename, "w") do f
        _save_dict_to_h5(f, my_data)
    end
end

### Reading from HDF5 files ###

function _load_dict_from_h5(group)
    data_dict = Dict()
    for key in keys(group)
        item = group[key]
        if isa(item, HDF5.Group)
            # If it's a group, recursively load it
            data_dict[key] = _load_dict_from_h5(item)
        else
            # Otherwise, it's a dataset
            data = read(item)
            # Julia automatically handles string decoding
            data_dict[key] = data
        end
    end
    return data_dict
end

function load_from_hdf5(filename)
    """
    Recursively load a nested dictionary of arrays from an HDF5 file.

    Parameters:
    filename (String): Name of the HDF5 file to read.
    
    Returns:
    Dict: A dictionary where keys are dataset names and values are arrays.
    """
    h5open(filename, "r") do f
        return _load_dict_from_h5(f)
    end
end

function load_key_from_hdf5(filename, key)
    """
    Load a specific key from an HDF5 file.

    Parameters:
    filename (String): Name of the HDF5 file to read.
    key (String): The key to load from the file.

    Returns:
    The data associated with the specified key.
    """
    h5open(filename, "r") do f
        if haskey(f, key)
            if isa(f[key], HDF5.Group)
                # If it's a group, recursively load it
                return _load_dict_from_h5(f[key])
            else
                return read(f[key])
            end
        else
            error("Key '$key' not found in the file '$filename'.")
        end
    end
end