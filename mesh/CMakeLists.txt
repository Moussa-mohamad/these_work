message(STATUS "________________________________CetSOLV : CMakeLists")
project (mesh)
message(STATUS "salut " ${CMAKE_MODULE_PATH})
find_package (Cython REQUIRED)
find_package (PythonLibs REQUIRED)


set (CYTHON_INCLUDE_MODULES numpy)

find_program (_PYTHON_BIN python NAMES python_d PATHS ${PYTHON_ROOT_DIR} NO_DEFAULT_PATH)

if (_PYTHON_BIN)
  foreach (module ${CYTHON_INCLUDE_MODULES})
  execute_process (COMMAND ${_PYTHON_BIN} -c "import numpy; exit (numpy.get_include())" OUTPUT_VARIABLE ${module}_INC_DIR ERROR_VARIABLE ${module}_INC_DIR)
  if (${module}_INC_DIR)
    include_directories (${${module}_INC_DIR})
  endif()
  endforeach()
endif()
unset (_PYTHON_BIN CACHE)

file (GLOB headers   *.hxx)
file (GLOB sources   *.cxx)
file (GLOB cyheaders *.pxd)
file (GLOB cysources *.pyx)

#WRAP_CYTHON_FILES (${cysources} FALSE cywrappers ${MESHGEMS_ROOT_DIR})
WRAP_CYTHON_FILES (${cysources} ON cywrappers "")

source_group ("Header Files"  FILES ${headers})
source_group ("Source Files"  FILES ${sources})
source_group ("Cython Files"  FILES ${cyheaders} ${cysources})
source_group ("Wrapper Files" FILES ${cywrappers})

link_directories (${INTEL_COMPILER_DIR}) # necessaire pour la commande suivante, raison non identifiee
python_add_module (mesh ${cywrappers} ${cyheaders} ${cysources} ${headers} ${sources})
set_target_properties (mesh PROPERTIES DEBUG_POSTFIX "_d")
set_target_properties (mesh PROPERTIES LINK_FLAGS "/NODEFAULTLIB:\"LIBCMT.lib\"")
target_link_libraries (mesh  
                                ${PYTHON_LIBRARY}
								"C:/Users/mmoussa/Desktop/gmsh-4.12.2-Windows64-sdk/lib/gmsh.lib"
								"C:/Users/mmoussa/Desktop/gmsh-4.12.2-Windows64-sdk/lib"
                                )
                                
message(STATUS '${PYTHON_LIBRARY}')
target_include_directories(mesh PUBLIC
                                          ${PYTHON_INCLUDE_DIR}
)
set_target_properties (mesh PROPERTIES COMPILE_FLAGS "/openmp")
set_property(TARGET mesh PROPERTY FOLDER "mesh")
install (TARGETS mesh LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})


