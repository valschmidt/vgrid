
#find_package( Boost COMPONENTS system filesystem REQUIRED )
#find_package(GDAL)
#include_directories(${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} "libMOOSGeodesy/include/MOOS/libMOOSGeodesy/MOOSGeodesy.h" ".")
#include_directories(${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} ".")

SET(SRC
    main.cpp
)

ADD_EXECUTABLE(vgrid_test ${SRC})

TARGET_LINK_LIBRARIES( vgrid_test
	vgrid
   ${Boost_SYSTEM_LIBRARY}
   ${Boost_FILESYSTEM_LIBRARY}
)

