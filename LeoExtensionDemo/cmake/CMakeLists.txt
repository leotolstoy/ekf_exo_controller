cmake_minimum_required(VERSION 2.8.1 FATAL_ERROR)

# # Setup Boilerplate Modification Variables
set(FormalProjectName AttitudeEstimator)
set(TestExecutable testAttitudeEstimator)
set(LibrarySourcesList 
 ../src/AttitudeEstimator.cpp 
 ../src/LowPassPositionFilter.cpp 
 ../src/HumePoseKalmanFilter.cpp 
 ../src/MovingPoseEstimator.cpp
 ../src/IMUReceiver.cpp)
set(TestSourcesList 
../test/TestAttitudeEstimator.cpp 
../test/TestAttitudeEstimator.hpp
../test/LiveTest.hpp
../test/LiveTest.cpp
../test/TestWithVelocity.hpp
../test/TestWithVelocity.cpp
../test/TestNoEigen.hpp
../test/TestNoEigen.cpp)
set(PublicIncludesList 
../include/AttitudeEstimator.hpp 
../include/LowPassPositionFilter.hpp 
../include/HumePoseKalmanFilter.hpp 
../include/MovingPoseEstimator.hpp
../include/noEigenAttitudeEstimatorInterface.hpp
../include/IMUReceiver.hpp
../include/IMUMessage.hpp)
  # TODO: separate library for eigen-less use?


# # Begin Project
project(${FormalProjectName})

# # Set Version
set (${FormalProjectName}_VERSION_MAJOR 0)
set (${FormalProjectName}_VERSION_MINOR 1)
set(${FormalProjectName}_PATCH_VERSION 0)
set(${FormalProjectName}_VERSION
  "${${FormalProjectName}_VERSION_MAJOR}.${${FormalProjectName}_VERSION_MINOR}.${${FormalProjectName}_PATCH_VERSION}")

# # Check for the GrayUtils library home environment variable
if(DEFINED ENV{GRAY_UTILS_HOME})
	message("gray utils home is defined $ENV{GRAY_UTILS_HOME}")
	set(CMAKE_MODULE_PATH "$ENV{GRAY_UTILS_HOME}/cmake" ${CMAKE_MODULE_PATH})
	set(GRAY_UTILS_CMAKE_FOLDER "$ENV{GRAY_UTILS_HOME}/cmake")
else()
	message(FATAL_ERROR "gray utils home doesn't exist. Please specify its location.")
endif()

# # Check for the GrayUtils library data environment variable
if(DEFINED ENV{GRAY_UTILS_DATA})
  message("gray utils data is defined $ENV{GRAY_UTILS_DATA}")
  set(CMAKE_HINTS_LOCATION "$ENV{GRAY_UTILS_DATA}/cmake")
else()
  message(FATAL_ERROR "gray utils home doesn't exist. Please specify its location.")
endif()

# # Set Installation Directories
set(INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib )
set(INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin )
set(INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/${FormalProjectName})
set(INSTALL_CMAKE_DIR ${CMAKE_HINTS_LOCATION} )
SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
 
find_package(Eigen3 REQUIRED)
find_package(GUThreads HINTS ${CMAKE_HINTS_LOCATION} REQUIRED CONFIGURE)
find_package(GUTiming HINTS ${CMAKE_HINTS_LOCATION} REQUIRED CONFIGURE)
find_package(PhaseSpaceUDPReceiver HINTS ${CMAKE_HINTS_LOCATION} REQUIRED CONFIGURE)
find_package(VectorGeometry HINTS ${CMAKE_HINTS_LOCATION} REQUIRED CONFIGURE)

# # Include the standard directorys for grayutils projects
include_directories(
  ${${FormalProjectName}_SOURCE_DIR}/../src 
  ${${FormalProjectName}_SOURCE_DIR}/../include  
  ${${FormalProjectName}_SOURCE_DIR}/../test
  ${GUThreads_INCLUDE_DIRS}
  ${GUTiming_INCLUDE_DIRS}
  ${PhaseSpaceUDPReceiver_INCLUDE_DIRS}
  ${VectorGeometry_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  )
message("EIGEN3_INCLUDE_DIR is ${EIGEN3_INCLUDE_DIR}")

# # Activate C++11 mode 
list( APPEND CMAKE_CXX_FLAGS "-std=c++11  ${CMAKE_CXX_FLAGS}")

# # Most projects build a library
add_library(${FormalProjectName} SHARED ${LibrarySourcesList} ${PublicIncludesList} )
target_link_libraries(${FormalProjectName} 
  ${GUThreads_LIBRARIES}
  ${GUTiming_LIBRARIES} 
  ${PhaseSpaceUDPReceiver_LIBRARIES}
  ${VectorGeometry_LIBRARIES}
  -lstdc++
  )


# # Here we have a test execuatable as well
add_executable(${TestExecutable} ${TestSourcesList} ${PublicIncludesList} )
target_link_libraries(${TestExecutable} 
  ${FormalProjectName} 
  -lstdc++
  )

# # It is important to associate the public includes with that library
set_target_properties(${FormalProjectName} PROPERTIES 
  OUTPUT_NAME "${FormalProjectName}"
	PUBLIC_HEADER "${PublicIncludesList}" 
	SOVERSION "${${FormalProjectName}_VERSION}")

# # Make a list of the targets you want to install
set(FullListOfTargets 
  ${FormalProjectName} 
  ${TestExecutable} 
  )

# #  Install the targets you have added in one go.
install (TARGETS ${FullListOfTargets} EXPORT ${FormalProjectName}Targets
	RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
	ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT archive
	LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
	PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev)

# # Add all targets to the build-tree export set
export(TARGETS ${FullListOfTargets}
  FILE "${CMAKE_BINARY_DIR}/${FormalProjectName}Targets.cmake")

# # Export the package, likely generating the <Proj>Targets-noconfig.cmake file
export(PACKAGE ${FormalProjectName})

# # For name space separation, preface configure replacements with CONF_
set(CONF_INCLUDE_DIRS "${INSTALL_INCLUDE_DIR}")
set(CONF_PROJ_NAME "${FormalProjectName}")
set(CONF_PROJ_LIBS ${FormalProjectName})
set(CONF_PROJ_BINS "")

# # Replace the @CONF_*@ symbols and rename the generic file.
configure_file(${GRAY_UTILS_CMAKE_FOLDER}/GenericConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${FormalProjectName}Config.cmake" @ONLY)

# # Use a generic config version template to handle version comparison logic
set(CONF_VERSION "${${FormalProjectName}_VERSION}")
configure_file(${GRAY_UTILS_CMAKE_FOLDER}/GenericConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/${FormalProjectName}ConfigVersion.cmake" @ONLY)

# # Install the two configured cmake files to the install cmake dir (grayutils/cmake)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${FormalProjectName}Config.cmake"
  "${PROJECT_BINARY_DIR}/${FormalProjectName}ConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
install(EXPORT ${FormalProjectName}Targets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)

install(CODE "execute_process(COMMAND sudo python setup.py install WORKING_DIRECTORY ../python_ext/)")

