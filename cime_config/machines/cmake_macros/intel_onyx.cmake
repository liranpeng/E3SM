string(APPEND FFLAGS " -fimf-use-svml")
if (NOT DEBUG)
  string(APPEND FFLAGS " -qno-opt-dynamic-align")
endif()
string(APPEND SLIBS " -lpthread")
