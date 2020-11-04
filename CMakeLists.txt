cmake_minimum_required(VERSION 3.13)
project(projet_ips1)
set(CMAKE_CXX_STANDARD 11)
set(COMPILE_OPTIONS PRIVATE -Wall -Wextra) #)-Werror -pedantic -ansi -Wshadow -Wdouble-promotion -Wundef -fno-common -Wconversion -Wunused-parameter)

#Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})


#MAIN
add_executable(main src/Hermite.cpp src/main.cpp src/SolverSchrodinger.cpp src/Derivator.cpp src/OrthogonalityChecker.cpp src/Saver.cpp)
target_link_libraries(main ${ARMADILLO_LIBRARIES})
target_compile_options(main ${COMPILE_OPTIONS})

#TESTS
add_executable(tests src/Hermite.cpp src/OrthogonalityChecker.cpp src/SolverSchrodinger.cpp src/Derivator.cpp
        tests/testsHermite.cpp tests/testsSolverSchrodinger.cpp tests/testsOrthogonalityChecker.cpp tests/testsDerivator.cpp)
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