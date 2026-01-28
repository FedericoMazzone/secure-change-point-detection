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
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_STATIC "Set to ON to include static versions of the library" OFF)
option(USE_SYSTEM_OPENFHE "Use OpenFHE already installed on the system" OFF)

set(OPENFHE_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/openfhe-development-1.1.2"
    CACHE PATH "Bundled OpenFHE source directory")

# Directory di build di OpenFHE (quando bundled)
set(OPENFHE_BUNDLED_BUILD_DIR
    "${CMAKE_BINARY_DIR}/_deps/openfhe-build"
    CACHE PATH "Bundled OpenFHE build directory")

# ------------------------------------------------------------
# OpenFHE: o di sistema (find_package) o bundled (add_subdirectory)
# ------------------------------------------------------------
if(USE_SYSTEM_OPENFHE)

    find_package(OpenFHE CONFIG REQUIRED)

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

    set(CMAKE_CXX_FLAGS ${OpenFHE_CXX_FLAGS})

    include_directories(${OPENMP_INCLUDES})
    include_directories(${OpenFHE_INCLUDE})
    include_directories(${OpenFHE_INCLUDE}/third-party/include)
    include_directories(${OpenFHE_INCLUDE}/core)
    include_directories(${OpenFHE_INCLUDE}/pke)
    include_directories(${OpenFHE_INCLUDE}/binfhe)

    link_directories(${OpenFHE_LIBDIR})
    link_directories(${OPENMP_LIBRARIES})

    if(BUILD_STATIC)
        set(CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
        link_libraries(${OpenFHE_STATIC_LIBRARIES})
    else()
        set(CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS})
        link_libraries(${OpenFHE_SHARED_LIBRARIES})
    endif()

else()

    message(STATUS "Building bundled OpenFHE from ${OPENFHE_SOURCE_DIR}")

    # Make OpenFHE self-contained (even if we never run "install")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_deps/openfhe-install" CACHE PATH "" FORCE)

    set(BUILD_SHARED ON CACHE BOOL "Build shared OpenFHE libraries" FORCE)
    set(BUILD_UNITTESTS OFF CACHE BOOL "Disable OpenFHE unit tests" FORCE)
    set(BUILD_EXAMPLES OFF CACHE BOOL "Disable OpenFHE examples" FORCE)
    set(BUILD_BENCHMARKS OFF CACHE BOOL "Disable OpenFHE benchmarks" FORCE)

    add_subdirectory("${OPENFHE_SOURCE_DIR}" "${OPENFHE_BUNDLED_BUILD_DIR}")


    # Include pragmatici (struttura standard OpenFHE)
    include_directories("${OPENFHE_SOURCE_DIR}/src")
    include_directories("${OPENFHE_SOURCE_DIR}/src/core")
    include_directories("${OPENFHE_SOURCE_DIR}/src/pke")
    include_directories("${OPENFHE_SOURCE_DIR}/src/binfhe")
    include_directories("${OPENFHE_SOURCE_DIR}/third-party/include")
    include_directories("${OPENFHE_SOURCE_DIR}/third-party/cereal/include")
    include_directories("${OPENFHE_SOURCE_DIR}/src/pke/include")
    include_directories("${OPENFHE_SOURCE_DIR}/src/core/include")
    include_directories("${OPENFHE_SOURCE_DIR}/src/binfhe/include")
    include_directories("${OPENFHE_BUNDLED_BUILD_DIR}/src/core")
    include_directories("${OPENFHE_BUNDLED_BUILD_DIR}/src/core/include")
    include_directories("${OPENFHE_BUNDLED_BUILD_DIR}/third-party/include")
    include_directories("${OPENFHE_BUNDLED_BUILD_DIR}/third-party/cereal/include")

    # Dove finiscono di solito le .so buildate (a seconda della piattaforma)
    set(OPENFHE_BUNDLED_LIBDIR1 "${OPENFHE_BUNDLED_BUILD_DIR}/lib")
    set(OPENFHE_BUNDLED_LIBDIR2 "${OPENFHE_BUNDLED_BUILD_DIR}/lib64")

endif()

# ------------------------------------------------------------
# Il tuo eseguibile
# ------------------------------------------------------------
set(SOURCES
    src/utils-basics.cpp
    src/utils-eval.cpp
    src/utils-matrices.cpp
    src/utils-ptxt.cpp
    src/cpd.cpp
)

add_executable(cpd ${SOURCES})

# Se OpenFHE è bundled, linka ai target OpenFHE e imposta RPATH
if(NOT USE_SYSTEM_OPENFHE)

    # Target OpenFHE tipici (se la tua versione usa nomi diversi, l'errore lo dirà chiaramente)
    target_link_libraries(cpd PRIVATE OPENFHEpke OPENFHEcore OPENFHEbinfhe)

    # RPATH: evita LD_LIBRARY_PATH quando si usa OpenFHE bundled
    set_target_properties(cpd PROPERTIES
        BUILD_RPATH   "${OPENFHE_BUNDLED_LIBDIR1};${OPENFHE_BUNDLED_LIBDIR2}"
        INSTALL_RPATH "${OPENFHE_BUNDLED_LIBDIR1};${OPENFHE_BUNDLED_LIBDIR2}"
    )

endif()
