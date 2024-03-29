from __future__ import annotations
import xarray as xr
import numpy as np
from scipy.spatial.transform import Rotation as R
from .constants import *

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
            if not rotation.single: 
                # If 'rotation' is a stack of rotations, apply each rotation sequentially
                for single_rotation in rotation:
                    vector = single_rotation.apply(vector)
                return vector
            else:
                # If 'rotation' represents a single rotation, apply it directly
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


    def rotate(self, rotvec=None, pole=None, skip_vars=['space','Ip']):
        """
        Rotates the coordinate system such that the z-axis is aligned with an input pole.  The new pole is defined by the input vector. 
        This will change all variables in the Dataset that have the "space" dimension, except for those passed to the skip_vars parameter.
        
        Parameters
        ----------
        ds : SwiftestDataset
            Dataset containing the vector quantity
        rotvec: (N,3) or (3,) float array
            Rotation vector
        pole : (3) float array
            New pole vector
        skip_vars : list of str, optional
            List of variable names to skip. The default is ['space','Ip'].
            
        Returns
        -------
        ds : SwiftestDataset
            Dataset with the new pole vector applied to all variables with the "space" dimension
            
        Notes
        -----
        You can pass either rotvec or pole, but not both. If both, or none, are passed, the function will raise an exception. 
        """
        if rotvec is not None and pole is not None:
            raise ValueError("You can only pass either rotvec or pole, but not both")
        if rotvec is None and pole is None:
            raise ValueError("You must pass either rotvec or pole")
         
        if 'space' not in self.dims:
            raise ValueError("Dataset must have a 'space' dimension")    
        
        # Verify that the new pole is a 3-element array
        if pole is not None:
            if len(pole) != 3:
                raise ValueError("Pole vector must be a 3-element array")
    
            # Normalize the new pole vector to ensure it is a unit vector
            pole_mag = np.linalg.norm(pole)
            unit_pole = pole / pole_mag
        
            # Define the original and target vectors
            target_vector = np.array([0, 0, 1])  # Rotate so that the z-axis is aligned with the new pole
            original_vector = unit_pole.reshape(1, 3)  
        
            # Use align_vectors to get the rotation that aligns the z-axis with Mars_rot
            rotvec, _ = R.align_vectors(target_vector, original_vector)
        elif rotvec is not None:
            rotvec = np.asarray(rotvec)
            if (rotvec.shape[-1]) != 3:
                raise ValueError("Rotation vector must be a 3-element array")
            
            rotvec = R.from_rotvec(rotvec) 

        # Loop through each variable in the dataset and apply the rotation if 'space' dimension is present
        for var in self.variables:
            if 'space' in self[var].dims and var not in skip_vars:
                self[var] = self[var].rotate(rotvec)
                
        return self
    
        
    def el2xv(self, GMcb: xr.DataArray | float | None = None) -> SwiftestDataset:
        """
        Converts a Dataset's orbital elements to Cartesian state vectors. The DataArray must have the appropriate dimensions for orbital elements.
        
        Parameters
        ----------
        GMcb : xr.DataArray or float
            Gravitational parameter of the central body
        
        Returns
        -------
        SwiftestDataset
            Dataset containing the computed state vectors (position 'rh' and velocity 'vh').
        """
        from .core import el2xv  # Assuming el2xv is implemented in the .core module

        if 'space' not in self.dims:
            raise ValueError("Dataset must have a 'space' dimension")
        
        required_vars = ['a', 'e', 'inc', 'capom', 'omega', 'capm']
        for var in required_vars:
            if var not in self.variables:
                raise ValueError(f"Dataset must have '{var}' variables")
        
        # Identify the index dimension
        if 'id' in self.dims:
            index_dim = 'id'
        elif 'name' in self.dims:
            index_dim = 'name'
        else:
            raise ValueError("Dataset must have an 'id' or 'name' dimension")
        
        if GMcb is None:
            if 'Gmass' not in self:
                raise ValueError("Dataset must have a 'Gmass' variable for the central body")
            if 'particle_type' in self.variables:
                GMcb = self['Gmass'].where(self['particle_type'] == CB_TYPE_NAME, drop=True)
            else:
                GMcb = self['Gmass'].where(self['id'] == 0, drop=True)
            if GMcb.size != 1:
                raise ValueError("Dataset must have a single central body") 
        
        if isinstance(GMcb, xr.DataArray):
            if 'id' in GMcb.dims:
                GMcb = GMcb.isel(id=0)
            elif 'name' in GMcb.dims:
                GMcb = GMcb.isel(name=0)
                
        if isinstance(GMcb, xr.DataArray):
            if 'id' in GMcb.dims:
                GMcb = GMcb.isel(id=0)
            elif 'name' in GMcb.dims:
                GMcb = GMcb.isel(name=0)
        else:
            GMcb = xr.DataArray(data = GMcb)
            
        for dim in self.dims:
            if dim not in GMcb.dims and dim not in ['space','l','m','sign']:
                GMcb = GMcb.expand_dims(dim={dim: self[dim]}) 
                
        if 'Gmass' in self:
            mu = xr.where(self['Gmass'] > 0.0, GMcb + self['Gmass'], GMcb)
        else:
            mu = GMcb
        
        # Prepare the orbital elements for the function call
        a = self['a']
        e = self['e']
        inc = self['inc']
        capom = self['capom']
        omega = self['omega']
        capm = self['capm']

        # Use apply_ufunc to convert orbital elements back to state vectors
        rh, vh = xr.apply_ufunc(
            el2xv,  # Function to apply
            mu, a, e, inc, capom, omega, capm,  # Inputs
            input_core_dims=[[index_dim], [index_dim], [index_dim], [index_dim], [index_dim], [index_dim], [index_dim]],  # Core dimensions for each input
            output_core_dims=[[index_dim, 'space'], [index_dim, 'space']],  # Core dimensions for outputs (position and velocity vectors)
            vectorize=True,  # Automatically vectorize over non-core dimensions
            dask="parallelized",  # Enable parallelized computation for Dask arrays, if applicable
            output_dtypes=[np.float64, np.float64]  # Expected data types for outputs
        )

        # Create a new Dataset with the state vectors
        new_vars = {'rh': rh, 'vh': vh}
        dataset = xr.Dataset(new_vars)
        if "name" in dataset.variables:
            dataset = dataset.drop_vars("name")        
        dsnew = xr.merge([dataset, self], compat="override")
        
        return SwiftestDataset(dsnew)
     
                
    def xv2el(self, GMcb: xr.DataArray | float | None = None) -> SwiftestDataset:
        """
        Converts A Dataset's Cartesian state vectors to orbital elements. The DataArray must have the "space" dimension. 
        
        Parameters
        ----------
        GMcb : xr.DataArray or float 
            Gravitational parameter of the central body
        
        Returns
        -------
        SwiftestDataset
            Dataset containing the computed orbital elements.
        """
        from .core import xv2el
        if 'space' not in self.dims:
            raise ValueError("Dataset must have a 'space' dimension")
        if 'rh' not in self.variables or 'vh' not in self.variables:
            raise ValueError("Dataset must have 'rh' and 'vh' variables")
       
        if 'id' in self.dims:
            index_dim = 'id'
        elif 'name' in self.dims:
            index_dim = 'name'
        else:
            raise ValueError("Dataset must have an 'id' or 'name' dimension")
        
        if GMcb is None:
            if 'Gmass' not in self:
                raise ValueError("Dataset must have a 'Gmass' variable for the central body")
            if 'particle_type' in self.variables:
                GMcb = self['Gmass'].where(self['particle_type'] == CB_TYPE_NAME, drop=True)
            else:
                GMcb = self['Gmass'].where(self['id'] == 0, drop=True)
            if GMcb.size != 1:
                raise ValueError("Dataset must have a single central body") 
        
        if isinstance(GMcb, xr.DataArray):
            if 'id' in GMcb.dims:
                GMcb = GMcb.isel(id=0)
            elif 'name' in GMcb.dims:
                GMcb = GMcb.isel(name=0)
        else:
            GMcb = xr.DataArray(data = GMcb)
            
        for dim in self.dims:
            if dim not in GMcb.dims and dim not in ['space','l','m','sign']:
                GMcb = GMcb.expand_dims(dim={dim: self[dim]}) 
                
        if 'Gmass' in self:
            mu = xr.where(self['Gmass'] > 0.0, GMcb + self['Gmass'], GMcb)
        else:
            mu = GMcb
        
        result = xr.apply_ufunc(
            xv2el,  
            mu,  
            self['rh'],  
            self['vh'],  
            input_core_dims=[[index_dim], [index_dim, 'space'], [index_dim, 'space']],  
            output_core_dims=[[index_dim]] * 11,  
            vectorize=True,  
            dask="parallelized",  
            output_dtypes=[np.float64] * 11  
        )
        
        varnames = ['a', 'e', 'inc', 'capom', 'omega', 'capm', 'varpi', 'lam', 'f', 'cape', 'capf']
        dataset = xr.Dataset({var: result[i] for i, var in enumerate(varnames)})
        if "name" in dataset.variables:
            dataset = dataset.drop_vars("name")        
        dsnew = xr.merge([dataset, self], compat="override")
        
        return SwiftestDataset(dsnew)

