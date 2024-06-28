# swiftest.add_body(**kwargs)
| Key Word Name   | Key Word Description                                                                                                                    | Options                        |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------|--------------------------------|
| ```name```      | Name(s) of bodies.                                                                                                                      | string or array-like of strings
| ```id```        | Unique identification value(s) of bodies.                                                                                               | float or array-like of floats
| ```a```         | Semi-major axis value(s) of bodies. Only used if  ```init_cond_format``` is set to ```EL```.                                            | float or array-like of floats
| ```e```         | Eccentricity value(s) of bodies. Only used if  ```init_cond_format``` is set to ```EL```.                                               | float or array-like of floats
| ```inc```       | Inclination value(s) of bodies in degrees. Only used if  ```init_cond_format``` is set to ```EL```.                                                | float or array-like of floats
| ```capom```     | Longitude of the ascending node value(s) of bodies in degrees. Only used if  ```init_cond_format``` is set to ```EL```.                            | float or array-like of floats
| ```omega```     | Argument of pericenter value(s) of bodies in degrees. Only used if  ```init_cond_format``` is set to ```EL```.                                     | float or array-like of floats
| ```capm```      | Mean anomaly value(s) of bodies in degrees. Only used if  ```init_cond_format``` is set to ```EL```.                                               | float or array-like of floats
| ```rh```        | Position vector(s) of bodies. Only used if  ```init_cond_format``` is set to ```XV```.                                                  | (n,3) array-like of floats
| ```vh```        | Velocity vector(s) of bodies. Only used if  ```init_cond_format``` is set to ```XV```.                                                  | (n,3) array-like of floats
| ```mass```      | Mass value(s) of bodies. Only for massive bodies. Only  ```mass``` **OR** ```Gmass``` may be set.                                       | float or array-like of floats
| ```Gmass```     | Gravitational mass value(s) of bodies. Only for massive bodies. Only  ```mass``` **OR** ```Gmass``` may be set.                         | float or array-like of floats
| ```radius```    | Radius value(s) of bodies. Only for massive bodies.                                                                                     | float or array-like of floats
| ```rhill```     | Hill Radius value(s) of bodies. Only for massive bodies.                                                                                | float or array-like of floats
| ```rot```       | Rotation rate vector(s) of bodies in degrees/TU. Only for massive bodies. Only used if ```rotation``` is set to ```True```.                           | (n,3) array-like of floats
| ```Ip```        | Principal axes moments of inertia vector(s) of bodies. Only for massive bodies. Only used if ```rotation``` is set to ```True```.       | (n,3) array-like of floats
| ```J2```        | The unitless value of the spherical harmonic term equal to J2*R^2 where R is the radius of the central body.                                                                                                     | float or array-like of floats
| ```J4```        | The unitless value of the spherical harmonic term equal to J4*R^4 where R is the radius of the central body.                                                                                                     | float or array-like of floats
