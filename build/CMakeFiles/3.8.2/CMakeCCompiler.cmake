set(CMAKE_C_COMPILER "/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpicc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "Intel")
set(CMAKE_C_COMPILER_VERSION "18.0.3.20180410")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_static_assert;c_variadic_macros;c_std_11")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_static_assert;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_SIMULATE_VERSION "")

set(CMAKE_AR "/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/ar")
set(CMAKE_RANLIB "/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/ranlib")
set(CMAKE_LINKER "/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/ld")
set(CMAKE_COMPILER_IS_GNUCC )
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "mpi;imf;svml;irng;m;ipgo;decimal;cilkrts;stdc++;irc;svml;c;irc_s;dl;c")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/dealii/8.5.0/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/gsl/2.5/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/p4est/1.1/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/ifort/2018.3.222/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64;/cvmfs/restricted.computecanada.ca/easybuild/software/2017/Core/ifort/2018.3.222/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64;/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/icc/2018.3.222/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64;/cvmfs/restricted.computecanada.ca/easybuild/software/2017/Core/icc/2018.3.222/compilers_and_libraries_2018.3.222/linux/tbb/lib/intel64/gcc4.4;/cvmfs/restricted.computecanada.ca/easybuild/software/2017/Core/icc/2018.3.222/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64;/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/gcc-7.3.0/lib64;/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/gcc-7.3.0/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/imkl/2018.3.222/mkl/lib/intel64;/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/imkl/2018.3.222/lib/intel64;/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib;/cvmfs/restricted.computecanada.ca/easybuild/software/2017/Core/icc/2018.3.222/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64_lin;/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib64;/cvmfs/soft.computecanada.ca/nix/store/rrwlh3bahkdwnbjvzm0nkq3504v451yl-gfortran-7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0;/cvmfs/soft.computecanada.ca/nix/store/rrwlh3bahkdwnbjvzm0nkq3504v451yl-gfortran-7.3.0/lib64;/cvmfs/soft.computecanada.ca/nix/store/rrwlh3bahkdwnbjvzm0nkq3504v451yl-gfortran-7.3.0/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
