cmake_minimum_required(VERSION 3.15)

### To use gcc/g++ on a Macintosh, you must set the Compilers
### here, not inside the project
##if(APPLE)
##       set(CMAKE_C_COMPILER "/usr/local/bin/gcc-7")
##       set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-7")
##endif()
### TODO: for now, we use CLang for Mac
###
### In order to create OpenFHE's static libraries you should enable
### the BUILD_STATIC option. For that, you run "cmake .. -DBUILD_STATIC=ON".
### After having your link completed you will find static libs
### with the suffix "_static" in ./build/libs/.
### Examples: OPENFHEpke_static.a, OPENFHEcore_static.a, etc.
### After you run "make install" in your build directory, you can build your custom application.
### If you need your application to be linked statically, then run "cmake .. -DBUILD_STATIC=ON"

project(demo CXX)
set(CMAKE_CXX_STANDARD 17)
option(BUILD_STATIC "Set to ON to include static versions of the library" OFF)
option(USE_SYSTEM_OPENFHE "Use OpenFHE already installed on the system" OFF)

set(OPENFHE_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/openfhe-development-1.1.2"
    CACHE PATH "Bundled OpenFHE source directory")

set(OPENFHE_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}/_deps/openfhe/install"
    CACHE PATH "Local OpenFHE install prefix")

include(ExternalProject)

if(NOT USE_SYSTEM_OPENFHE)
    
    find_package(OpenFHE CONFIG QUIET)

    if(NOT OpenFHE_FOUND)
        
        message(STATUS "OpenFHE not found. Building bundled OpenFHE...")

        ExternalProject_Add(openfhe_ep
            SOURCE_DIR  "${OPENFHE_SOURCE_DIR}"
            PREFIX      "${CMAKE_BINARY_DIR}/_deps/openfhe"
            INSTALL_DIR "${OPENFHE_INSTALL_PREFIX}"

            CMAKE_ARGS
                -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                -DBUILD_SHARED=ON
                -DBUILD_UNITTESTS=OFF
                -DBUILD_EXAMPLES=OFF

            BUILD_COMMAND   ${CMAKE_COMMAND} --build <BINARY_DIR> -j
            INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR>
        )

        list(PREPEND CMAKE_PREFIX_PATH "${OPENFHE_INSTALL_PREFIX}")

        find_package(OpenFHE CONFIG REQUIRED)

    endif()

else()

    find_package(OpenFHE CONFIG REQUIRED)

endif()

if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif ()

set( CMAKE_CXX_FLAGS ${OpenFHE_CXX_FLAGS} )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${OpenFHE_INCLUDE} )
include_directories( ${OpenFHE_INCLUDE}/third-party/include )
include_directories( ${OpenFHE_INCLUDE}/core )
include_directories( ${OpenFHE_INCLUDE}/pke )
include_directories( ${OpenFHE_INCLUDE}/binfhe )
### add directories for other OpenFHE modules as needed for your project

link_directories( ${OpenFHE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} )
endif()

### ADD YOUR EXECUTABLE(s) HERE
### add_executable( EXECUTABLE-NAME SOURCES )
###
### EXAMPLE:
### add_executable( test demo-simple-example.cpp )

set( SOURCES
        src/utils-basics.cpp
        src/utils-eval.cpp
        src/utils-matrices.cpp
        src/utils-ptxt.cpp
        src/cpd.cpp
)

add_executable( cpd ${SOURCES} )

if(NOT USE_SYSTEM_OPENFHE)
    set_target_properties(cpd PROPERTIES
        BUILD_RPATH   "${OPENFHE_INSTALL_PREFIX}/lib"
        INSTALL_RPATH "${OPENFHE_INSTALL_PREFIX}/lib"
    )
endif()

if(TARGET openfhe_ep)
    add_dependencies(cpd openfhe_ep)
endif()