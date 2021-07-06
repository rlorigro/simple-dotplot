######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.10)
project(simple_dotplot VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")

# Compiler options.
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS )
#set(ASAN_OPTIONS=check_initialization_order=1)
#set(ASAN_OPTIONS=detect_leaks=1)

set(CMAKE_CXX_FLAGS "-fsanitize=address -fexceptions")  # leak checking and address sanitization
#set(CMAKE_CXX_FLAGS "-fexceptions")

#add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
add_definitions(-O3 -Wall)              # Much optimization

# Definitions needed to eliminate runtime dependency
# on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")


#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/Color.cpp
        src/FastaReader.cpp
        src/FastaWriter.cpp
        src/FastqIterator.cpp
        src/Path
        )

project(simple_dotplot)
add_library(simple_dotplot STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(simple_dotplot PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(simple_dotplot PROPERTIES DEFINE_SYMBOL "")

##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- mummer4 ----------------- #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# Download or update library as an external project
ExternalProject_Add(project_mummer
        URL https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
#        PREFIX ${CMAKE_SOURCE_DIR}/external/mummer
#        DOWNLOAD_COMMAND ""
#        UPDATE_COMMAND ""
#        STAMP_DIR ${CMAKE_SOURCE_DIR}/external/mummer/stamp
        BUILD_IN_SOURCE TRUE
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/mummer/
        CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/external/mummer/configure --prefix=${CMAKE_SOURCE_DIR}/external/mummer/
        INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/mummer/
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_mummer INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(mummer4 STATIC IMPORTED)
set_property(TARGET mummer4
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libumdmummer.a)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(mummer4
        project_mummer)

add_dependencies(simple_dotplot mummer4)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/mummer/include/)

# Add runtime path for main library so it can access external library
#set_property(TARGET GetBlunted PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/ogdf")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################


# -------- TESTS --------

set(TESTS
        test_mummer
        test_cli
        )

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            simple_dotplot
            mummer4
            )

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        dotplot
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            simple_dotplot
            mummer4
            )

endforeach()

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        simple_dotplot
        DESTINATION ${BINARY_INSTALL_DIR})
