cmake_minimum_required(VERSION 3.10)

# Set the source and header directories
set(SOURCES_DIR "${CMAKE_SOURCE_DIR}/../Sources")
set(HEADER_DIR "${CMAKE_SOURCE_DIR}/../Header Files")


# Set the build directory for final output
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../build")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/CMakeFiles")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/CMakeFiles")
set(CMAKE_CACHEFILE_DIR "${CMAKE_SOURCE_DIR}/CMakeFiles")

# Add executable target
add_executable(Triangulation
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
)

# Include header files
target_include_directories(Triangulation PUBLIC 
  ${HEADER_DIR}
  ${PREVIOUS_PROJECT_HEADER_DIR}
)

# Find CGAL
find_package(CGAL QUIET COMPONENTS Qt5)

# Find Boost
find_package(Boost REQUIRED COMPONENTS system filesystem)

# Find Qt5
find_package(Qt5 REQUIRED COMPONENTS Core Gui)

# Set up the compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 -g")
set(CMAKE_BUILD_TYPE Debug)

# Link against CGAL, Boost, and Qt5
target_link_libraries(Triangulation
  CGAL::CGAL
  CGAL::CGAL_Qt5
  Boost::system
  Boost::filesystem
  Qt5::Core
  Qt5::Gui
)

# Add the definitions for CGAL
add_definitions(-DCGAL_USE_BASIC_VIEWER)
