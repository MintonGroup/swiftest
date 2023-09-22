set(CMAKE_Fortran_COMPILER "/opt/homebrew/bin/gfortran-12")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "12.3.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_Fortran_SIMULATE_VERSION "")


set(CMAKE_Fortran_COMPILER_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX13.sdk/usr")
set(CMAKE_COMPILER_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX13.sdk/usr")

set(CMAKE_AR "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/opt/homebrew/bin/gcc-ar-12")
set(CMAKE_RANLIB "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib")
set(CMAKE_TAPI "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/tapi")
set(CMAKE_Fortran_COMPILER_RANLIB "/opt/homebrew/bin/gcc-ranlib-12")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95;f03;F03;f08;F08)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
set(CMAKE_Fortran_LINKER_DEPFILE_SUPPORTED FALSE)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()


set(CMAKE_Fortran_SYSROOT_FLAG "-isysroot")
set(CMAKE_Fortran_OSX_DEPLOYMENT_TARGET_FLAG "-mmacosx-version-min=")

set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/opt/homebrew/Cellar/gcc@12/12.3.0/lib/gcc/12/gcc/aarch64-apple-darwin22/12/finclude;/usr/local/include;/opt/homebrew/include;/opt/homebrew/Cellar/gcc@12/12.3.0/lib/gcc/12/gcc/aarch64-apple-darwin22/12/include;/opt/homebrew/Cellar/gcc@12/12.3.0/lib/gcc/12/gcc/aarch64-apple-darwin22/12/include-fixed;/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;emutls_w;gcc;quadmath;emutls_w;gcc;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/usr/local/lib;/opt/homebrew/lib;/opt/homebrew/Cellar/gcc@12/12.3.0/lib/gcc/12/gcc/aarch64-apple-darwin22/12;/opt/homebrew/Cellar/gcc@12/12.3.0/lib/gcc/12/gcc;/opt/homebrew/Cellar/gcc@12/12.3.0/lib/gcc/12")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
