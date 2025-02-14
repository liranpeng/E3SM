if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_SLASHPROC")
endif()
set(ESMF_LIBDIR "/projects/ccsm/esmf-6.3.0rp1/lib/libO/Linux.intel.64.openmpi.default")
if (MPILIB STREQUAL openmpi)
  set(MPI_PATH "$ENV{MPIHOME}")
endif()
if (MPILIB STREQUAL mpi-serial AND NOT compile_threaded)
  set(PFUNIT_PATH "/projects/ccsm/pfunit/3.2.9/mpi-serial")
endif()
set(PIO_FILESYSTEM_HINTS "lustre")
