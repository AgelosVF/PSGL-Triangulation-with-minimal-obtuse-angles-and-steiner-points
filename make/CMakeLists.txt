cmake_minimum_required(VERSION 3.10)
project(opt_triangulation)

# Set Release as default build type and disable dev warnings
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
    set(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION OFF CACHE BOOL "" FORCE)
    # Disable developer warnings
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON CACHE INTERNAL "")
    message(STATUS "Setting build type to Release as none was specified")
endif()

#old cmake policy for find_package
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

# Set the source and header directories
set(SOURCES_DIR "${CMAKE_SOURCE_DIR}/../Sources")
set(HEADER_DIR "${CMAKE_SOURCE_DIR}/../Header Files")


# Set the build directory for final output
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../build")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/Lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/Lib")
set(CMAKE_CACHEFILE_DIR "${CMAKE_SOURCE_DIR}")

# Add executable target
add_executable(opt_triangulation
  ${SOURCES_DIR}/main.cpp
  ${SOURCES_DIR}/boost_utils.cpp
  ${SOURCES_DIR}/triangulation_utils.cpp
  ${SOURCES_DIR}/SteinerPoints.cpp
  ${SOURCES_DIR}/outputfile.cpp
  ${SOURCES_DIR}/LocalSearch.cpp
  ${SOURCES_DIR}/SimulatedAnnealing.cpp
  ${SOURCES_DIR}/Previous_Project.cpp
  ${SOURCES_DIR}/AntCollony.cpp
  ${SOURCES_DIR}/BoundaryType.cpp
  ${SOURCES_DIR}/Random.cpp
)

# Include header files
target_include_directories(opt_triangulation PUBLIC 
  ${HEADER_DIR}
  ${PREVIOUS_PROJECT_HEADER_DIR}
)

# Find CGAL with better error handling
find_package(CGAL QUIET COMPONENTS Qt5)
if(NOT CGAL_FOUND)
    message(FATAL_ERROR "CGAL not found!")
endif()

# Find Boost with better error handling
find_package(Boost REQUIRED COMPONENTS system filesystem)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost not found!")
endif()

# Find Qt5 with better error handling
find_package(Qt5 REQUIRED COMPONENTS Core Gui)
if(NOT Qt5_FOUND)
    message(FATAL_ERROR "Qt5 not found!")
endif()

# Set up the compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 -g")


# Link against CGAL, Boost, and Qt5
target_link_libraries(opt_triangulation
  CGAL::CGAL
  CGAL::CGAL_Qt5
  Boost::system
  Boost::filesystem
  Qt5::Core
  Qt5::Gui
)

# Add the definitions for CGAL
add_definitions(-DCGAL_USE_BASIC_VIEWER)
