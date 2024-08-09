cmake_minimum_required (VERSION 3.20) 

project(panmanUtils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations") 
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

find_package(Protobuf REQUIRED CONFIG) 
include_directories(${Protobuf_INCLUDE_DIRS})
if(Protobuf_FOUND)
	message(STATUS "Found Protocol Buffers: ${Protobuf_VERSION}")
endif()

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

if(DEFINED Protobuf_PATH)
	add_library(panmanUtils STATIC
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
	
	add_library(panmanUtils STATIC
		${PANMAT_SRCS}
		${PROTO_SRCS}
		${PROTO_HDRS}
		)

endif()

TARGET_COMPILE_OPTIONS(panmanUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

TARGET_LINK_LIBRARIES(panmanUtils PRIVATE stdc++ JsonCpp::JsonCpp ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES})
target_include_directories(panmanUtils PUBLIC "${PROJECT_BINARY_DIR}")

# Add installation rules
install(TARGETS panmanUtils
        EXPORT panmanUtilsTargets
        ARCHIVE DESTINATION lib)

install(DIRECTORY src/
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")

# Export the targets to a script
install(EXPORT panmanUtilsTargets
        FILE panmanUtilsTargets.cmake
        NAMESPACE panmanUtils::
        DESTINATION lib/cmake/panmanUtils)

# Create a Config file for find_package
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/panmanUtilsConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_file(panmanUtilsConfig.cmake.in
               "${CMAKE_CURRENT_BINARY_DIR}/panmanUtilsConfig.cmake"
               @ONLY)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/panmanUtilsConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/panmanUtilsConfigVersion.cmake"
        DESTINATION lib/cmake/panmanUtils)

# Export from the build tree
export(EXPORT panmanUtilsTargets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/panmanUtilsTargets.cmake"
       NAMESPACE panmanUtils::)

export(PACKAGE panmanUtils)