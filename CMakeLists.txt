project(MPC)

cmake_minimum_required (VERSION 3.5)

# It's better to use -DCMAKE_BUILD_TYPE=Release
# for building optimized executable
if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++1z")

file(GLOB_RECURSE sources src/mpc/*.cpp src/mpc/*.h)

include_directories("${CMAKE_SOURCE_DIR}/src")

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        include_directories(/usr/local/opt/openssl/include)
        link_directories(/usr/local/opt/openssl/lib)
        link_directories(/usr/local/Cellar/libuv/1.11.0/lib)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

add_executable(mpc ${sources})

target_link_libraries(mpc ipopt z ssl uv uWS)

