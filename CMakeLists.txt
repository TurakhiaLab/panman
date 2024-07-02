cmake_minimum_required (VERSION 3.8) 

project(panmanUtils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -lboost_filesystem -lboost_system -lboost_program_options") 
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

if(DEFINED Protobuf_PATH) 
	find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 
else()
	find_package(Protobuf REQUIRED) 
endif()

# Print version messages
if(Protobuf_FOUND)
	message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
endif()
include_directories(${Protobuf_INCLUDE_DIRS})

# Include TBB
include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


# Include JSONCPP
include(${CMAKE_TOOLCHAIN_FILE})
find_package(jsoncpp CONFIG REQUIRED)


file(GLOB PANMAT_SRCS "src/main.cpp" "src/panmanUtils.cpp" "src/panmanUtils.hpp" "src/panman.hpp")

if(DEFINED Protobuf_PATH)
	add_executable(panmanUtils
		${PANMAT_SRCS}
		)

	protobuf_generate(
		LANGUAGE cpp
		TARGET panmanUtils 
		PROTOS panman.proto)


else()
	protobuf_generate_cpp(
		PROTO_SRCS PROTO_HDRS
		panman.proto)
	
	add_executable(panmanUtils
		${PANMAT_SRCS}
		${PROTO_SRCS}
		${PROTO_HDRS}
		)

endif()

TARGET_COMPILE_OPTIONS(panmanUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

TARGET_LINK_LIBRARIES(panmanUtils PRIVATE stdc++ JsonCpp::JsonCpp ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) #${Boost_LIBRARIES}  ) # OpenMP::OpenMP_CXX)
target_include_directories(panmanUtils PUBLIC "${PROJECT_BINARY_DIR}")
