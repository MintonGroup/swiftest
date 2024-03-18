import xarray as xr
import numpy as np
from scipy.spatial.transform import Rotation as R

class SwiftestDataArray(xr.DataArray):
    """
    N-dimensional ``xarray.DataArray``-like array. Inherits from ``xarray.DataArray`` and has its own set of methods and attributes 
    specific to the Swiftest project

    Parameters
    ----------------
    *args:
        Arguments for the ``xarray.DataArray`` class
    **kwargs:
        Keyword arguments for the ``xarray.DataArray`` class

    Notes
    -----
    See `xarray.DataArray <https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html>`__
    for further information about DataArrays.
    """
    
    __slots__ = ()

    def __init__(self, *args,**kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def _construct_direct(cls, *args, **kwargs):
        """Override to make the result a ``swiftest.SwiftestDataArray`` class."""
        return cls(xr.DataArray._construct_direct(*args, **kwargs))

    def to_dataset(self): 
        """Converts a ``SwiftestDataArray`` into a ``SwiftestDataset`` with a single data
        variable."""
        xrds = super().to_dataset()
        return SwiftestDataset(xrds)
    
    def magnitude(self, name: str | None = None): 
        """
        Computes the magnitude of a vector quantity. Note: The DataArray must have the "space" dimension. 
        
        Parameters
        ----------
        name : str, optional
            Name of the new DataArray. By default, the string "_mag" is appended to the original name.
        
        Returns
        -------
        mag : SwiftestDataArray
            DataArray containing the magnitude of the vector quantity 
        """
        dim = "space"
        ord = None
        
        if dim not in self.dims:
            raise ValueError(f"Dimension {dim} not found in DataArray")
        if name is None:
            name = self.name + "_mag"
        da = xr.apply_ufunc(
            np.linalg.norm, self.where(~np.isnan(self)), input_core_dims=[[dim]], kwargs={"ord": ord, "axis": -1}, dask="allowed"
        )
        
        da = da.rename(name)
        return SwiftestDataArray(da)


    def rotate(self, rotation):
        """
        Rotates a vector quantity using a rotation matrix. The DataArray must have the "space" dimension. 
        
        Parameters
        ----------
        rotation : (3) float array
            Rotation vector
        """
       
        if "space" not in self.dims:
            raise ValueError("DataArray must have a 'space' dimension")
        
        # Define a function to apply the rotation, which will be used with apply_ufunc
        def apply_rotation(vector, rotation):
            return rotation.apply(vector)        
        
        da = xr.apply_ufunc(
            apply_rotation,
            self,
            kwargs={'rotation': rotation},
            input_core_dims=[['space']],
            output_core_dims=[['space']],
            vectorize=True,
            dask='parallelized',
            output_dtypes=[self.dtype]
        )
        return SwiftestDataArray(da) 
    

class SwiftestDataset(xr.Dataset):
    """
    A ``xarray.Dataset``-like, multi-dimensional, in memory, array database.  Inherits from ``xarray.Dataset`` and has its own set of
    methods and attributes specific to the Swiftest project.

    Parameters
    ----------------
    *args:
        Arguments for the ``xarray.Dataset`` class
    **kwargs:
        Keyword arguments for the ``xarray.Dataset`` class

    Notes
    -----
    See `xarray.Dataset <https://docs.xarray.dev/en/stable/generated/xarray.Dataset.html>`__
    for further information about Datasets.
    """
    __slots__ = ()
    def __init__(
                 self,
                 *args,
                 **kwargs,
                ):
        super().__init__(*args, **kwargs)

    def __getitem__(self, key):
        """Override to make sure the result is an instance of
        ``swiftest.SwiftestDataArray`` or ``swiftest.SwiftestDataset``."""

        value = super().__getitem__(key)

        if isinstance(value, xr.DataArray):
            value = SwiftestDataArray(value)
        elif isinstance(value, xr.Dataset):
            value = SwiftestDataset(value)

        return value

    def _calculate_binary_op(self, *args, **kwargs):
        """Override to make the result a complete instance of ``swiftest.SwiftestDataset``."""
        ds = super()._calculate_binary_op(*args, **kwargs)
        if not isinstance(ds, SwiftestDataset):
            ds = SwiftestDataset(ds)
        ds = SwiftestDataset(ds)

        return ds

    def _construct_dataarray(self, name) -> SwiftestDataArray:
        """Override to make the result an instance of ``swiftest.SwiftestDataArray``."""
        xarr = super()._construct_dataarray(name)
        return SwiftestDataArray(xarr)

    @classmethod
    def _construct_direct(cls, *args, **kwargs):
        """Override to make the result an ``swiftest.SwiftestDataset`` class."""

        return cls(xr.Dataset._construct_direct(*args, **kwargs))


    def _replace(self, *args, **kwargs):
        """Override to make the result a complete instance of ``swiftest.SwiftestDataset``."""
        ds = super()._replace(*args, **kwargs)

        if not isinstance(ds, SwiftestDataset):
            ds = SwiftestDataset(ds)

        return ds

    @classmethod
    def from_dataframe(cls, dataframe):
        """Override to make the result a ``swiftest.SwiftestDataset`` class."""

        return cls(
            {col: ("index", dataframe[col].values) for col in dataframe.columns},
            coords={"index": dataframe.index},
        )

    @classmethod
    def from_dict(cls, data, **kwargs):
        """Override to make the result a ``swiftest.SwiftestDataset`` class."""

        return cls(
            {key: ("index", val) for key, val in data.items()},
            coords={"index": range(len(next(iter(data.values()))))},
            **kwargs,
        )


    def rotate_to_pole(ds, new_pole, skip_vars=['space','Ip']):
        """
        Rotates the coordinate system such that the z-axis is aligned with an input pole.  The new pole is defined by the input vector. 
        This will change all variables in the Dataset that have the "space" dimension, except for those passed to the skip_vars parameter.
        
        Parameters
        ----------
        ds : SwiftestDataset
            Dataset containing the vector quantity
        new_pole : (3) float array
            New pole vector
        skip_vars : list of str, optional
            List of variable names to skip. The default is ['space','Ip'].
            
        Returns
        -------
        ds : SwiftestDataset
            Dataset with the new pole vector applied to all variables with the "space" dimension
        """
        
        if 'space' not in ds.dims:
            raise ValueError("Dataset must have a 'space' dimension")    
        
        # Verify that the new pole is a 3-element array
        if len(new_pole) != 3:
            print("New pole must be a 3-element array")
            return ds
    
        # Normalize the new pole vector to ensure it is a unit vector
        pole_mag = np.linalg.norm(new_pole)
        unit_pole = new_pole / pole_mag
        
        # Define the original and target vectors
        target_vector = np.array([0, 0, 1])  # Rotate so that the z-axis is aligned with the new pole
        original_vector = unit_pole.reshape(1, 3)  
        
        # Use align_vectors to get the rotation that aligns the z-axis with Mars_rot
        rotvec, _ = R.align_vectors(target_vector, original_vector)

        # Loop through each variable in the dataset and apply the rotation if 'space' dimension is present
        for var in ds.variables:
            if 'space' in ds[var].dims and var not in skip_vars:
                ds[var] = ds[var].rotate(rotvec)
                
        return ds
            
                