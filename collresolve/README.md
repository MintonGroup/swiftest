**collresolve** is a library designed to provide collision analysis and handling for N-body codes.

The library has interfaces for the following languages:
* C/C++
* Fortran
* Python

## Installation

### C/C++/Fortran

The installation of the C/C++/Fortran library using the standard `./configure`, `make` and `make install` commands.

The repository does not provide the `configure` script. If installing directly from the repository, the file must be created first. The easiest way to do this is executing `autoreconf --install` command which will create the `configure` script and its dependencies. This command requires the [GNU autoconf](http://www.gnu.org/software/autoconf) tool. In this situation, the commands to execute are:
```
autoreconf --install
./configure
make
make install
```

As for other packages, you may want to execute `./configure --help` to see which options are available and tune that command to your needs, e.g. by changing the location where the library will be installed using the `--prefix=PATH` option.

### Linking

If the library has been installed in a non-standard part, then the paths to the header files and library object files must be provided to the compiler and linker calls. Hereafter, we assume that the base location of the library is `PREFIX`, which is the value of the `--prefix` argument to `./configure` call.

The public interface of the library for the C/C++ languages is provided in the `collresolve.h` file. In case the path must be provided, then the argument `-IPREFIX/include` must be added to the compiler commands of files that make use of `collresolve.h`.

To use the library, the argument `-lcollresolve` is to be provided to the linking command. In case the library is in a non-standard directory, the path can be provided using `-LPREFIX/lib` argument.

The Fortran interface is tailored to use with the `mercury` package.

### Python

From a checkout of the repository, the Python module can be easily built and installed by the executing the following commands:
```
python setup.py build
python setup.py install
```

In case the library is to be installed for the current user only, then the `--user` argument can be provided to the `install` command.

This will put the Python module in a location where is it readily available. No further action is needed.

## Usage

To use most of the library, a configuration object must be created and set first. The configuration object contains essential parameters for the calls, such as the collision model to determine the outcome of collisions and the unit system in use. The object must not be accessed directly, but the functions `collresolve_conf_*` should be used instead to alter its state.

### Python

An example of usage of the library is provided in the `example.py` file.

### Mercury

To use the library with the `mercury`, the following modifications to the code of the latter are needed. When compiling the code, the additional flags described above to link the executable to the library are required.

In the `mce_coll` subroutine, near the end, in lieu of the part
```
c
c Do the collision (inelastic merger)
      call mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost,
     %         nsetup)
```
the code should be changed to something like
```
      if (opt(2).eq.2 .and. i.gt.1 .and. j.gt.1) then
        model = 4
        nres = 2

        regime = collresolve_resolve(model, m(i) / K2, m(j) / K2,
     %     rphys(i), rphys(j), xh(:,i), xh(:,j), vh(:,i), vh(:,j), nres,
     %     mres, rres, pres, vres)

        if (regime .lt. 0) then
          ! An error occurred
          ! Do not do anything.
        else if (mres(1) .lt. 1.d-3 / (1047.d0 * 317.8d0)) then
          stat(i) = -2
          stat(j) = -2
          xh(:,j) = -xh(:,j)
          vh(:,j) = -vh(:,j)
        else if (mres(2) .lt. 1.d-3 / (1047.d0 * 317.8d0)) then
          m(i) = mres(1) * K2
          rphys(i) = rres(1)
          xh(:,i) = pres(:,1)
          vh(:,i) = vres(:,1)
          stat(j) = -2
          xh(:,j) = -xh(:,j)
          vh(:,j) = -vh(:,j)
        else
          m(i) = mres(1) * K2
          m(j) = mres(2) * K2
          rphys(i) = rres(1)
          rphys(j) = rres(2)
          xh(:,i) = pres(:,1)
          xh(:,j) = pres(:,2)
          vh(:,i) = vres(:,1)
          vh(:,j) = vres(:,2)
        end if
      else
c
c Do the collision (inelastic merger)
        call mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost,
     %         nsetup)
      end if
```
with the following new variables at the beginning of the subroutine
```
      real*8 mres(3),rres(3),pres(3,3),vres(3,3)
      integer nres,model,regime,collresolve_resolve
```
There are a few items to be noted with the above code:
* The library is only called when the flag about using the "fragmentation" mode in the parameters file is enabled. This allows to easily perform comparison run with mercury's standard merging algorithm with having to re-compile the code.
* The value of the `model` variable should be adjusted to which model the library is to use to resolve the collisions. The possible values are given in the `collresolve_model` enummeration in `collresolve.h`. The value provided in this code snippet, `4`, tells the library to use Cambioni et al. (2019) model.
* It implements a minimum mass cutoff for the remnants, with a value to 1/1000 of an Earth mass (the factor 1047 * 317.8 being the conversion to solar mass that Mercury uses as the mass unit).

While the above codes sets the correct radii on return, these will be ignored in a standard version of the `mercury` package. In effect, this will assume a constant density for each body. In case the model 4 (Cambioni et al. 2019) is used, a further modification is needed so that the bodies have a consistent radius with the bodies that were used to generate the model. To achieve this, the code setting the physical radius in `mce_init` should be changed from `rphys(j)=hill(j)/a(j)*(temp/rho(j))**THIRD` to `rphys(j) = collresolve_radius(4, m(j) / K2)` while adding `real*8 collresolve_radius` in the definitions of that subroutine. This will make mercury use the library's mass-radius relation  for all bodies (both at the beginning of the simulation and after a collision), except for the central body. In effect, this makes the `d` parameter in the input file useless.

## License

The library is licensed under version 2.0 of the Apache License, see the `LICENSE` file for the full terms and conditions.

## Citations

If you use this library in a scientific work that lead to publication, we would like you to acknowledge the following article:
* Emsenhuber, A., Cambioni S., Asphaug, E., Gabriel, T. S. J., Schwartz, S. R., and Furfaro, R. (subm.). Realistic On-the-fly Outcomes of Planetary Collisions II: Bringing Machine Learning to N-body Simulations. **The Astrophysical Journal**.

If you use the `LS2012` model, you should also cite:
 * Leinhardt, Z. M. and Stewart, S. T. (2012). Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling Laws. **The Astrophysical Journal**, 745(1), 79. [doi:10.1088/0004-637X/745/1/79](https://doi.org/10.1088/0004-637X/745/1/79) [bib:2012ApJ...745...79L](https://ui.adsabs.harvard.edu/abs/2012ApJ...745...79L)

If you use the `SL2012` model, you should also cite the same publication as for `LS2012`, and:
 * Stewart, S. T. and Leinhardt, Z. M. (2012). Collisions between Gravity-dominated Bodies. II. The Diversity of Impact Outcomes during the End Stage of Planet Formation. **The Astrophysical Journal**, 751(1), 32. [doi:10.1088/0004-637X/751/1/32](https://doi.org/10.1088/0004-637X/751/1/32) [bib:2012ApJ...751...32S](https://ui.adsabs.harvard.edu/abs/2012ApJ...751...32S)
 * Genda, H., Kokubo, E., and Ida, S. (2012). Merging Criteria for Giant Impacts of Protoplanets. **The Astrophysical Journal**, 744(2), 137. [doi:10.1088/0004-637X/744/2/137](https://doi.org/10.1088/0004-637X/744/2/137) [bib:2012ApJ...744..137G](https://ui.adsabs.harvard.edu/abs/2012ApJ...744..137G)

If you use the `C2019` model, you should also cite:
 * Cambioni, S., Asphaug, E., Emsenhuber, A., Gabriel, T. S. J., Furfaro, R., and Schwartz, S. R. (2019). Realistic On-the-fly Outcomes of Planetary Collisions: Machine Learning Applied to Simulations of Giant Impacts. **The Astrophysical Journal**, 875(1), 40. doi:[10.3847/1538-4357/ab0e8a](https://doi.org/10.3847/1538-4357/ab0e8a) [bib:2019ApJ...875...40C](https://ui.adsabs.harvard.edu/abs/2019ApJ...875...40C)
