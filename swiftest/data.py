import xarray as xr
from typing import IO

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

     