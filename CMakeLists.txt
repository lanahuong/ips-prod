cmake_minimum_required(VERSION 3.13)
project(projet_ips1)
set(CMAKE_CXX_STANDARD 11)
set(COMPILE_OPTIONS PRIVATE -pedantic -ansi -Wshadow -Wdouble-promotion -Wundef -fno-common  -Wunused-parameter -Winit-self -Wuninitialized -Wmissing-declarations -Woverloaded-virtual -Wold-style-cast )
#-Wconversion

set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -mtune=native -fopenmp -Ofast -flto")

#Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})


#MAIN
add_executable(main
        src/Poly.cpp src/Poly.h
        src/main.cpp
        src/Basis.h
        src/Saver.cpp src/Saver.h
        src/constants.h
        src/Basis.cpp
        src/NuclearDensityCalculator.cpp src/NuclearDensityCalculator.h
        src/Chrono.hpp src/ThreadSafeAccumulator.hpp src/FactorisationHelper.hpp)
target_link_libraries(main ${ARMADILLO_LIBRARIES})
target_compile_options(main ${COMPILE_OPTIONS})

#TESTS
add_executable(tests src/Poly.cpp src/Basis.cpp src/Poly.h src/Basis.h
        src/Saver.cpp src/Saver.h src/constants.h
        src/NuclearDensityCalculator.cpp src/NuclearDensityCalculator.h
        tests/testsMandatory.cpp
        tests/testsNuclearDensityCalculator.cpp src/Chrono.hpp src/ThreadSafeAccumulator.hpp src/FactorisationHelper.hpp)
target_link_libraries(tests ${ARMADILLO_LIBRARIES})
target_compile_options(tests ${COMPILE_OPTIONS})
target_link_libraries(tests gtest_main)


# GOOGLE TEST
# from https://github.com/google/googletest/blob/master/googletest/README.md
# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if (result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if (result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif ()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif ()