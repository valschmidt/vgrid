#--------------------------------------------------------
# The CMakeLists.txt for:                       lib_vgrid
# Author(s):                                  Val Schmidt
#--------------------------------------------------------

# Boost and its components
find_package( Boost COMPONENTS system filesystem REQUIRED )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()

FILE(GLOB SRC *.cpp *.c)

SET(HEADERS
	vgrid.h
)

ADD_LIBRARY(vgrid ${SRC})

TARGET_LINK_LIBRARIES(vgrid
   ${Boost_SYSTEM_LIBRARY}
   ${Boost_FILESYSTEM_LIBRARY}
)

install(TARGETS vgrid DESTINATION lib)
install(FILES vgrid.h DESTINATION include)
