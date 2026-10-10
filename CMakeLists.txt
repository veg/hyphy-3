cmake_minimum_required(VERSION 3.20)
project(hyphy3 VERSION 3.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra")

# Header-only Eigen 3.4.0
include(FetchContent)
FetchContent_Declare(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_GetProperties(Eigen3)
if(NOT eigen3_POPULATED)
    FetchContent_Populate(Eigen3)
endif()

# Core interface library
add_library(hyphy3_core INTERFACE)
add_library(hyphy_core ALIAS hyphy3_core)

target_include_directories(hyphy3_core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${eigen3_SOURCE_DIR}
)

# OpenMP support
find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND AND EXISTS "/opt/homebrew/opt/libomp")
    set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
    find_package(OpenMP)
endif()
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP CXX found: ${OpenMP_CXX_VERSION}")
    target_link_libraries(hyphy3_core INTERFACE OpenMP::OpenMP_CXX)
endif()

# Test targets
enable_testing()
add_executable(hyphy_test tests/test_core.cpp)
target_link_libraries(hyphy_test PRIVATE hyphy3_core)
add_test(NAME CoreTests COMMAND hyphy_test)

add_executable(hyphy_fel_test tests/test_fel.cpp)
target_link_libraries(hyphy_fel_test PRIVATE hyphy3_core)
add_test(NAME FELTests COMMAND hyphy_fel_test)

add_executable(hyphy_fel_cd2_parity_test tests/test_fel_cd2_parity.cpp)
target_link_libraries(hyphy_fel_cd2_parity_test PRIVATE hyphy3_core)
add_test(NAME CD2ParityTest COMMAND hyphy_fel_cd2_parity_test)

add_executable(hyphy_fel_multi_dataset_test tests/test_fel_multi_dataset.cpp)
target_link_libraries(hyphy_fel_multi_dataset_test PRIVATE hyphy3_core)
add_test(NAME MultiDatasetParityTest COMMAND hyphy_fel_multi_dataset_test)

add_executable(hyphy_meme_test tests/test_meme.cpp)
target_link_libraries(hyphy_meme_test PRIVATE hyphy3_core)
add_test(NAME MEMEParityTest COMMAND hyphy_meme_test)

add_executable(hyphy_busted_test tests/test_busted.cpp)
target_link_libraries(hyphy_busted_test PRIVATE hyphy3_core)
add_test(NAME BUSTEDTest COMMAND hyphy_busted_test)

add_executable(hyphy_busted_mh_test tests/test_busted_mh.cpp)
target_link_libraries(hyphy_busted_mh_test PRIVATE hyphy3_core)
add_test(NAME BUSTEDMHTest COMMAND hyphy_busted_mh_test)

add_executable(hyphy_cpp_autograd_test tests/test_cpp_autograd.cpp)
target_link_libraries(hyphy_cpp_autograd_test PRIVATE hyphy3_core)
add_test(NAME CppAutogradTest COMMAND hyphy_cpp_autograd_test)

add_executable(hyphy_absrel_test tests/test_absrel.cpp)
target_link_libraries(hyphy_absrel_test PRIVATE hyphy3_core)
add_test(NAME ABSRELTest COMMAND hyphy_absrel_test)

add_executable(hyphy_relax_test tests/test_relax.cpp)
target_link_libraries(hyphy_relax_test PRIVATE hyphy3_core)
add_test(NAME RELAXTest COMMAND hyphy_relax_test)


# Unified HyPhy 3 CLI Driver (hyphy3 fel ..., hyphy3 meme ..., hyphy3 busted ..., hyphy3 absrel ..., hyphy3 relax ...)
add_executable(hyphy3
    src/apps/hyphy3.cpp
    src/apps/hyphy_fel.cpp
    src/apps/hyphy_meme.cpp
    src/apps/hyphy_busted.cpp
    src/apps/hyphy_absrel.cpp
    src/apps/hyphy_relax.cpp
)
target_compile_definitions(hyphy3 PRIVATE HYPHY3_COMBINED_DRIVER)
target_link_libraries(hyphy3 PRIVATE hyphy3_core)

# Standalone CLI Applications
add_executable(hyphy3-fel src/apps/hyphy_fel.cpp)
target_link_libraries(hyphy3-fel PRIVATE hyphy3_core)

add_executable(hyphy3-meme src/apps/hyphy_meme.cpp)
target_link_libraries(hyphy3-meme PRIVATE hyphy3_core)

add_executable(hyphy3-busted src/apps/hyphy_busted.cpp)
target_link_libraries(hyphy3-busted PRIVATE hyphy3_core)

add_executable(hyphy3-absrel src/apps/hyphy_absrel.cpp)
target_link_libraries(hyphy3-absrel PRIVATE hyphy3_core)

add_executable(hyphy3-relax src/apps/hyphy_relax.cpp)
target_link_libraries(hyphy3-relax PRIVATE hyphy3_core)

# Backward-compatible aliases
add_executable(hyphy-fel ALIAS hyphy3-fel)
add_executable(hyphy-meme ALIAS hyphy3-meme)
add_executable(hyphy-busted ALIAS hyphy3-busted)
add_executable(hyphy-absrel ALIAS hyphy3-absrel)
add_executable(hyphy-relax ALIAS hyphy3-relax)

# Python & nanobind bindings
find_package(Python 3.8 COMPONENTS Interpreter Development.Module)
if(Python_FOUND)
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_DIR
        RESULT_VARIABLE nanobind_STATUS
    )
    if(nanobind_STATUS EQUAL 0 AND EXISTS "${nanobind_DIR}")
        list(APPEND CMAKE_PREFIX_PATH "${nanobind_DIR}")
        find_package(nanobind CONFIG)
        if(nanobind_FOUND)
            message(STATUS "Building Python bindings with nanobind")
            nanobind_add_module(_hyphy3 NB_STATIC
                src/python/bindings.cpp
            )
            target_link_libraries(_hyphy3 PRIVATE hyphy3_core)
            set_target_properties(_hyphy3 PROPERTIES
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python"
            )

            add_test(NAME PythonAutogradTest
                COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/python/tests/test_python_autograd.py"
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            )
            set_tests_properties(PythonAutogradTest PROPERTIES
                ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python;KMP_DUPLICATE_LIB_OK=TRUE"
            )

            add_test(NAME JSONSchemaTest
                COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_json_schema.py"
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            )
            set_tests_properties(JSONSchemaTest PROPERTIES
                ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python;KMP_DUPLICATE_LIB_OK=TRUE"
            )
        endif()
    endif()
endif()

