if (NOT DEBUG)
  string(APPEND CFLAGS " -O2 -kind=byte")
endif()
string(APPEND CPPDEFS " -DLINUX")
if (NOT DEBUG)
  string(APPEND FFLAGS " -kind=byte")
endif()
if (DEBUG)
  string(APPEND FFLAGS " -O0 -v")
endif()
set(PIO_FILESYSTEM_HINTS "lustre")
string(APPEND SLIBS " -lpmi")
