cmake_minimum_required(VERSION 3.18)

project(panmanUtils)

# ------------------ Versioning ------------------
set(panmanUtils_VERSION "?.?.?")
execute_process(COMMAND head -n1 ${CMAKE_SOURCE_DIR}/version.txt OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
if("${VERSION_PATCH}" MATCHES "[0-9]+")
    set(panmanUtils_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()
configure_file(src/version.hpp.in ${CMAKE_SOURCE_DIR}/src/version.hpp)

# ------------------ Build Flags ------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function -Wno-deprecated-declarations -w")


set(CMAKE_INCLUDE_CURRENT_DIR ON)

# ------------------ Dependencies ------------------

# Protobuf
if(DEFINED Protobuf_PATH) 
    find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH})
else()
    find_package(Protobuf REQUIRED)
endif()
message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
include_directories(${Protobuf_INCLUDE_DIRS})

# CapnProto
find_program(CAPNP_EXECUTABLE capnp REQUIRED)
find_program(CAPNPC_CXX_EXECUTABLE capnpc-c++ REQUIRED)
find_package(CapnProto REQUIRED)
message(STATUS "Using Capn Proto ${CapnProto_VERSION}")
include_directories(${CAPNP_INCLUDE_DIRS})
add_definitions(${CAPNP_DEFINITIONS})

# Boost
find_package(Boost REQUIRED COMPONENTS program_options iostreams filesystem date_time)
include_directories(${Boost_INCLUDE_DIRS})

# JSONCPP
if(CMAKE_TOOLCHAIN_FILE)
    include(${CMAKE_TOOLCHAIN_FILE})
endif()
find_package(jsoncpp CONFIG REQUIRED)

# Accept an installed classic TBB package, retaining the legacy source build.
if(EXISTS "${TBB_DIR}/cmake/TBBBuild.cmake")
    include(${TBB_DIR}/cmake/TBBBuild.cmake)
    tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
    find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)
else()
    find_package(TBB CONFIG REQUIRED)
endif()

# ------------------ Source Files ------------------
file(GLOB PANMAT_SRCS
    "src/panmanUtils.cpp"
    "src/panman.cpp"
    "src/panmanUtils.hpp"
    "src/panman.hpp"
)

# ------------------ Executable ------------------
add_executable(panmanUtils ${PANMAT_SRCS} src/constructionBackend.cpp)
target_include_directories(panmanUtils PRIVATE "${CMAKE_SOURCE_DIR}/src")
option(PANMAN_ENABLE_CUDA "Build optional CUDA construction inference" ON)
if(PANMAN_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        target_sources(panmanUtils PRIVATE gpu/constructionInference.cu)
        if(NOT CMAKE_CUDA_ARCHITECTURES)
            if(CMAKE_VERSION VERSION_LESS "3.23")
                message(FATAL_ERROR "Set CMAKE_CUDA_ARCHITECTURES explicitly or use CMake >= 3.23 for CUDA")
            endif()
            set(CMAKE_CUDA_ARCHITECTURES all-major)
        endif()
        set_target_properties(panmanUtils PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON
                              CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
        target_compile_definitions(panmanUtils PRIVATE PANMAN_HAS_CUDA=1)
        message(STATUS "Construction backend: CPU and CUDA")
    else()
        message(STATUS "CUDA compiler unavailable; construction auto uses CPU")
    endif()
endif()

# CapnProto and Protobuf Codegen
if(DEFINED CapnProto_PATH)
    capnp_generate(LANGUAGE cpp TARGET panmanUtils PROTOS panman.capnp)
    protobuf_generate(LANGUAGE cpp TARGET panmanUtils PROTOS usher.proto panman.proto)
else()
    capnp_generate_cpp(CAPNP_SRCS CAPNP_HDRS panman.capnp)
    protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS usher.proto panman.proto)
    target_sources(panmanUtils PRIVATE ${CAPNP_SRCS} ${CAPNP_HDRS} ${PROTO_SRCS} ${PROTO_HDRS})
endif()

target_compile_options(panmanUtils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

# Only link Abseil on macOS
if(APPLE AND Protobuf_VERSION VERSION_GREATER_EQUAL "22")
    find_package(absl CONFIG REQUIRED)
    target_link_libraries(panmanUtils PRIVATE
        absl::log
        absl::strings
        absl::base
    )
	# target_link_libraries(panmanUtils PRIVATE TBB::tbb)
endif()

# ------------------ Linking ------------------
target_link_libraries(panmanUtils PRIVATE
    JsonCpp::JsonCpp
    ${Boost_LIBRARIES}
    TBB::tbb
    ${TBB_IMPORTED_TARGETS}
    ${CAPNP_LIBRARIES}
    ${Protobuf_LIBRARIES}
)

target_include_directories(panmanUtils PUBLIC "${PROJECT_BINARY_DIR}")

# Small deterministic unit tests plus file-format integration tests. CUDA parity
# is exercised by the integration suite when a device is actually available.
include(CTest)
if(BUILD_TESTING)
    add_executable(panman_layout_tests tests/test_layout.cpp)
    target_include_directories(panman_layout_tests PRIVATE "${CMAKE_SOURCE_DIR}/src")
    add_test(NAME layout_and_inference COMMAND panman_layout_tests)
    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME maf_export COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_maf_export.py" "$<TARGET_FILE:panmanUtils>")
        add_test(NAME maf_multichrom COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_maf_multichrom.py" "$<TARGET_FILE:panmanUtils>")
        add_test(NAME construction_compatibility COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_construction.py" "$<TARGET_FILE:panmanUtils>" --capnp "${CAPNP_EXECUTABLE}")
    endif()
endif()
