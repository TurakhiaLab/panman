cmake_minimum_required (VERSION 3.8) 

project(panmanUtils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -lboost_filesystem -lboost_system -lboost_program_options") 
set(CMAKE_INCLUDE_CURRENT_DIR ON) 
set(BUILD_SHARED_LIBS ON)

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
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm")
	set(ARCH "arm64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "armv4t")
	set(ARCH "arm64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86")
	set(ARCH "x86")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
	set(ARCH "x86_64")
endif()

tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1 arch=${ARCH})
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview HINTS ${TBB_DIR})

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


# Include JSONCPP
include(${CMAKE_TOOLCHAIN_FILE})
find_package(jsoncpp CONFIG REQUIRED)


file(GLOB PANMAT_SRCS "src/panmanUtils.cpp" "src/panman.cpp" "src/panmanUtils.hpp" "src/panman.hpp")


	protobuf_generate_cpp(
		PROTO_SRCS PROTO_HDRS
		panman.proto)
	
	add_executable(panmanUtils
		${PANMAT_SRCS}
		${PROTO_SRCS}
		${PROTO_HDRS}
		)

target_include_directories(panmanUtils PUBLIC ${TBB_DIR}/include)

TARGET_COMPILE_OPTIONS(panmanUtils PUBLIC -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

TARGET_LINK_LIBRARIES(panmanUtils PUBLIC JsonCpp::JsonCpp ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} stdc++) #${Boost_LIBRARIES}  ) # OpenMP::OpenMP_CXX)
target_include_directories(panmanUtils PUBLIC "${PROJECT_BINARY_DIR}")


add_library(panman STATIC ${PANMAT_SRCS} ${PROTO_SRCS} ${PROTO_HDRS})
target_include_directories(panman PUBLIC "${PROJECT_BINARY_DIR}")
target_include_directories(panman PUBLIC ${TBB_DIR}/include)

target_link_libraries(panman PUBLIC JsonCpp::JsonCpp ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} stdc++)
