CMAKE_MINIMUM_REQUIRED(VERSION 3.6)
PROJECT( vgrid )

set(CMAKE_CXX_STANDARD 11)

# Add the libraries in the current directory to the include path
#============================================================================
FILE(GLOB LOCAL_LIBRARY_DIRS ./lib_*)
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/lib_vgrid")
LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/lib_vgrid")
find_package( Boost COMPONENTS system filesystem REQUIRED )

#============================================================================
# List the subdirectories to build...
#============================================================================
ADD_SUBDIRECTORY(lib_vgrid)
ADD_SUBDIRECTORY(vgrid)

