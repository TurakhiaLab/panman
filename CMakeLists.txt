cmake_minimum_required (VERSION 3.8) 

project(panmanUtils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -lboost_filesystem -lboost_system -lboost_program_options") 
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

# if(DEFINED Protobuf_PATH) 
# 	find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 
# else()
# 	find_package(Protobuf REQUIRED) 
# endif()

# Print version messages
# if(Protobuf_FOUND)
# 	message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
# endif()
# include_directories(${Protobuf_INCLUDE_DIRS})

# Adding capn proto
if(DEFINED CapnProto_PATH) 
	find_package(CapnProto REQUIRED HINTS ${CapnProto_PATH}) 
else()
	find_package(CapnProto REQUIRED) 
endif()

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
file(GLOB CAP_SRCS "src/capnprototest.cpp")

if(DEFINED CapnProto_PATH)
	add_executable(panmanUtils
		${PANMAT_SRCS}
		)

	# protobuf_generate(
	# 	LANGUAGE cpp
	# 	TARGET panmanUtils 
	# 	PROTOS panman.proto)
	
	capnp_generate(
		LANGUAGE cpp
		TARGET panmanUtils  
		PROTOS panman.capnp)

	add_executable(capnp ${CAP_SRCS})
	capnp_generate(
		LANGUAGE cpp
		TARGET capnp
		PROTOS panman.capnp
	)


else()
	# protobuf_generate_cpp(
	# 	PROTO_SRCS PROTO_HDRS
	# 	panman.proto)
	
	capnp_generate_cpp(
		CAPNP_SRCS CAPNP_HDRS 
		panman.capnp)
	
	add_executable(panmanUtils
		${PANMAT_SRCS}
		${CAPNP_SRCS}
		${CAPNP_HDRS}
		)
	
	add_executable(capnp
	${CAP_SRCS}
	${CAPNP_SRCS}
	${CAPNP_HDRS}
	)

endif()

TARGET_COMPILE_OPTIONS(panmanUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

TARGET_LINK_LIBRARIES(panmanUtils PRIVATE stdc++ JsonCpp::JsonCpp ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${CAPNP_LIBRARIES} ) #${Protobuf_LIBRARIES} ${Boost_LIBRARIES}  ) # OpenMP::OpenMP_CXX)
target_include_directories(panmanUtils PUBLIC "${PROJECT_BINARY_DIR}")

TARGET_COMPILE_OPTIONS(capnp PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

TARGET_LINK_LIBRARIES(capnp PRIVATE stdc++ JsonCpp::JsonCpp ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${CAPNP_LIBRARIES} ) #${Protobuf_LIBRARIES} ${Boost_LIBRARIES}  ) # OpenMP::OpenMP_CXX)
target_include_directories(capnp PUBLIC "${PROJECT_BINARY_DIR}")
