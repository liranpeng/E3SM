string(APPEND CONFIG_ARGS " --host=cray")
set(USE_CUDA "TRUE")
string(APPEND CPPDEFS " -DGPU")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")
string(APPEND CUDA_FLAGS " -ccbin CC -O2 -arch sm_80 --use_fast_math")
set(CXX_LINKER "FORTRAN")
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")
