cmake_minimum_required (VERSION 3.10) 

project(panmanUtils)

# set version
set(panmanUtils_VERSION "?.?.?")
execute_process(COMMAND head -n1 ${CMAKE_SOURCE_DIR}/version.txt OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
if ("${VERSION_PATCH}" MATCHES "[0-9]+")
    set(panmanUtils_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()
configure_file(src/version.hpp.in ${CMAKE_SOURCE_DIR}/src/version.hpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -w") 
if(NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f")
endif()
set(CMAKE_INCLUDE_CURRENT_DIR ON) 



# add google proto
if(DEFINED Protobuf_PATH) 
	find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 
else()
	find_package(Protobuf REQUIRED) 
endif()
if(Protobuf_FOUND)
	message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
endif()
include_directories(${Protobuf_INCLUDE_DIRS})


# add capn proto
find_program(CAPNP_EXECUTABLE capnp REQUIRED)
find_program(CAPNPC_CXX_EXECUTABLE capnpc-c++ REQUIRED)
find_package(CapnProto REQUIRED)

if(CapnProto_FOUND)
	message(STATUS "Using Capn Proto ${CapnProto_VERSION}")
endif()
include_directories(${CAPNP_INCLUDE_DIRS})
add_definitions(${CAPNP_DEFINITIONS})

# Include TBB
include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


# Include JSONCPP
include(${CMAKE_TOOLCHAIN_FILE})
find_package(jsoncpp CONFIG REQUIRED)


file(GLOB PANMAT_SRCS "src/panmanUtils.cpp" "src/panman.cpp" "src/panmanUtils.hpp" "src/panman.hpp")

if(DEFINED CapnProto_PATH)
	add_executable(panmanUtils
		${PANMAT_SRCS}
		)

	capnp_generate(
		LANGUAGE cpp
		TARGET panmanUtils  
		PROTOS panman.capnp)
	
	protouf_generate(
		LANGUAGE cpp
		TARGET panmanUtils  
		PROTOS usher.proto panman.proto)

else()
	capnp_generate_cpp(
		CAPNP_SRCS CAPNP_HDRS 
		panman.capnp)
	
	protobuf_generate_cpp(
		PROTO_SRCS PROTO_HDRS
		usher.proto panman.proto)
	
	add_executable(panmanUtils
		${PANMAT_SRCS}
		${CAPNP_SRCS}
		${CAPNP_HDRS}
		${PROTO_SRCS}
		${PROTO_HDRS}
		)
	
endif()

TARGET_COMPILE_OPTIONS(panmanUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

TARGET_LINK_LIBRARIES(panmanUtils PRIVATE stdc++ JsonCpp::JsonCpp ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${CAPNP_LIBRARIES} ${Protobuf_LIBRARIES}) #${Protobuf_LIBRARIES} ${Boost_LIBRARIES}  ) # OpenMP::OpenMP_CXX)
target_include_directories(panmanUtils PUBLIC "${PROJECT_BINARY_DIR}")
