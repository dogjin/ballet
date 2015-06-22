cmake_minimum_required ( VERSION 2.6.4 )

project ( Ballet )

# package related variables
set (PACKAGE "ballet")
set (PACKAGE_NAME "Ballet Signal Processing & Simulation Platform")
set (PACKAGE_VERSION "1.0.1")
set (PACKAGE_VERSION_MAJOR 1)
set (PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set (PACKAGE_TARNAME ${PACKAGE})

add_subdirectory ( ballet )
