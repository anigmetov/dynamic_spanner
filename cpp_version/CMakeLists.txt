project                     (spanner_wasserstein)
cmake_minimum_required      (VERSION 3.5.1)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include (GenerateExportHeader)
find_package(Armadillo REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost COMPONENTS system filesystem REQUIRED)

# Default to Release

if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)

include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                             ${CMAKE_CURRENT_SOURCE_DIR}/include/cover_tree
                             ${CMAKE_CURRENT_SOURCE_DIR}/include/wasserstein
                             ${ARMADILLO_INCLUDE_DIR}
                             SYSTEM ${Boost_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 14)

if(NOT WIN32)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${OpenMP_CXX_FLAGS} -Wall -pedantic -Wextra -pedantic-errors")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -D_GLIBCXX_DEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -O2 -g -ggdb")
endif(NOT WIN32)

file(GLOB WS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/wasserstein/*.h ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB CT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/cover_tree/*.h)
file(GLOB DS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

file(GLOB SRC_FILES  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

find_package (Threads)
set (libraries ${libraries}     ${CMAKE_THREAD_LIBS_INIT} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${ARMADILLO_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

add_executable(spanner_wasserstein ${SRC_FILES} ${WS_HEADERS} ${CT_HEADERS} ${DS_HEADERS})
target_link_libraries(spanner_wasserstein PUBLIC ${libraries} )
