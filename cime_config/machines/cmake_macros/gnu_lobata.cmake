set(SUPPORTS_CXX "TRUE")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_VPRINTF -DHAVE_TIMES -DHAVE_GETTIMEOFDAY -DHAVE_BACKTRACE")
endif()
string(APPEND SLIBS " -lcurl")
if (DEBUG)
  string(APPEND FFLAGS " -g -fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow")
else()
  string(APPEND FFLAGS " -fno-unsafe-math-optimizations")
endif()
if (MPILIB STREQUAL mpi-serial)
  set(SCC "gcc")
  set(SFC "gfortran")
endif()
string(APPEND CXX_LIBS " -lstdc++")
