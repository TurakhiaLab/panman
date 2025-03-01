cmake_minimum_required (VERSION 3.8) 

project(panmanUtils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -lboost_filesystem -lboost_system -lboost_program_options -w") 
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

find_package(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


# Include JSONCPP
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
