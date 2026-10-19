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

# HTML results viewer (viz/fel-viewer), embedded into the binaries at build time so hyphy3 can
# write self-contained reports with no Node, Python or companion files (cmake/EmbedViewer.cmake)
set(HYPHY3_VIEWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/viz/fel-viewer)
set(HYPHY3_OHMIC_VIEWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/viz/ohmic-viewer)
set(HYPHY3_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(HYPHY3_VIEWER_HEADER ${HYPHY3_GENERATED_DIR}/hyphy/report/viewer_assets.hpp)
set(HYPHY3_VIEWER_STAMP ${HYPHY3_GENERATED_DIR}/viewer_assets.stamp)
file(GLOB_RECURSE HYPHY3_VIEWER_SOURCES CONFIGURE_DEPENDS
    ${HYPHY3_VIEWER_DIR}/src/*.js
    ${HYPHY3_VIEWER_DIR}/src/*.css
    ${HYPHY3_OHMIC_VIEWER_DIR}/src/*.js
    ${HYPHY3_OHMIC_VIEWER_DIR}/src/*.css
)
add_custom_command(
    OUTPUT ${HYPHY3_VIEWER_STAMP}
    BYPRODUCTS ${HYPHY3_VIEWER_HEADER}
    COMMAND ${CMAKE_COMMAND}
        -DVIEWER_DIR=${HYPHY3_VIEWER_DIR}
        -DOHMIC_VIEWER_DIR=${HYPHY3_OHMIC_VIEWER_DIR}
        -DOUTPUT=${HYPHY3_VIEWER_HEADER}
        -DSTAMP=${HYPHY3_VIEWER_STAMP}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmbedViewer.cmake
    DEPENDS ${HYPHY3_VIEWER_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmbedViewer.cmake
    COMMENT "Embedding HTML results viewer"
    VERBATIM
)
add_custom_target(hyphy3_viewer_assets DEPENDS ${HYPHY3_VIEWER_STAMP})

add_library(hyphy3_report INTERFACE)
target_include_directories(hyphy3_report INTERFACE ${HYPHY3_GENERATED_DIR})
target_link_libraries(hyphy3_report INTERFACE hyphy3_core)
add_dependencies(hyphy3_report hyphy3_viewer_assets)

# Test targets
enable_testing()
# Lets tests locate benchmarks/ and tests/data/ regardless of working directory
add_compile_definitions(HYPHY3_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
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

add_executable(hyphy_contrast_fel_test tests/test_contrast_fel.cpp)
target_link_libraries(hyphy_contrast_fel_test PRIVATE hyphy3_core)
add_test(NAME ContrastFELTest COMMAND hyphy_contrast_fel_test)

add_executable(hyphy_fubar_test tests/test_fubar.cpp)
target_link_libraries(hyphy_fubar_test PRIVATE hyphy3_core)
add_test(NAME FUBARTest COMMAND hyphy_fubar_test)

add_executable(hyphy_prime_test tests/test_prime.cpp)
target_link_libraries(hyphy_prime_test PRIVATE hyphy3_core)
add_test(NAME PRIMETest COMMAND hyphy_prime_test)

add_executable(hyphy_gard_test tests/test_gard.cpp)
target_link_libraries(hyphy_gard_test PRIVATE hyphy3_core)
add_test(NAME GARDTest COMMAND hyphy_gard_test)

add_executable(hyphy_slac_test tests/test_slac.cpp)
target_link_libraries(hyphy_slac_test PRIVATE hyphy3_core)
add_test(NAME SLACTest COMMAND hyphy_slac_test)

add_executable(hyphy_dfel_test tests/test_dfel.cpp)
target_link_libraries(hyphy_dfel_test PRIVATE hyphy3_core)
add_test(NAME DirectionalFELTest COMMAND hyphy_dfel_test)

add_executable(hyphy_ohmic_test tests/test_ohmic.cpp)
target_link_libraries(hyphy_ohmic_test PRIVATE hyphy3_core)
add_test(NAME OhmicTest COMMAND hyphy_ohmic_test)

add_executable(hyphy_report_test tests/test_report.cpp)
target_link_libraries(hyphy_report_test PRIVATE hyphy3_report)
add_test(NAME ReportTest COMMAND hyphy_report_test)

add_executable(hyphy_expcm_test tests/test_expcm.cpp)
target_link_libraries(hyphy_expcm_test PRIVATE hyphy3_core)
add_test(NAME ExpcmTest COMMAND hyphy_expcm_test)

add_executable(hyphy_expcm_optimizer_test tests/test_expcm_optimizer.cpp src/analyses/expcm_optimizer.cpp)
target_link_libraries(hyphy_expcm_optimizer_test PRIVATE hyphy3_core)
add_test(NAME ExpcmOptimizerTest COMMAND hyphy_expcm_optimizer_test)

add_executable(hyphy_optimizer_quick_test tests/test_optimizer_quick.cpp src/analyses/expcm_optimizer.cpp)
target_link_libraries(hyphy_optimizer_quick_test PRIVATE hyphy3_core)
add_test(NAME OptimizerQuickTest COMMAND hyphy_optimizer_quick_test)

add_executable(hyphy_phase4_persite_test tests/test_phase4_persite.cpp)
target_link_libraries(hyphy_phase4_persite_test PRIVATE hyphy3_core)
add_test(NAME Phase4PersiteTest COMMAND hyphy_phase4_persite_test)

add_executable(hyphy_phase4_optimizer_test tests/test_phase4_optimizer.cpp src/analyses/persite_optimizer.cpp)
target_link_libraries(hyphy_phase4_optimizer_test PRIVATE hyphy3_core)
add_test(NAME Phase4OptimizerTest COMMAND hyphy_phase4_optimizer_test)

add_executable(hyphy_phase4_validation_test tests/test_phase4_validation.cpp src/analyses/persite_optimizer.cpp)
target_link_libraries(hyphy_phase4_validation_test PRIVATE hyphy3_core)
add_test(NAME Phase4ValidationTest COMMAND hyphy_phase4_validation_test)

add_test(NAME HTMLReportEndToEndTest
    COMMAND ${CMAKE_COMMAND}
        -DHYPHY3=$<TARGET_FILE:hyphy3>
        -DDATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/data
        -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/html_report_test
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_html_report.cmake
)


# Unified HyPhy 3 CLI Driver (hyphy3 fel ..., hyphy3 meme ..., hyphy3 busted ..., hyphy3 absrel ..., hyphy3 relax ..., hyphy3 fubar ..., hyphy3 prime ..., hyphy3 gard ..., hyphy3 slac ..., hyphy3 dfel ..., hyphy3 ohmic ...)
add_executable(hyphy3
    src/apps/hyphy3.cpp
    src/apps/hyphy_fel.cpp
    src/apps/hyphy_contrast_fel.cpp
    src/apps/hyphy_meme.cpp
    src/apps/hyphy_busted.cpp
    src/apps/hyphy_absrel.cpp
    src/apps/hyphy_relax.cpp
    src/apps/hyphy_fubar.cpp
    src/apps/hyphy_prime.cpp
    src/apps/hyphy_gard.cpp
    src/apps/hyphy_slac.cpp
    src/apps/hyphy_dfel.cpp
    src/apps/hyphy_ohmic.cpp
    src/apps/hyphy_report.cpp
)
target_compile_definitions(hyphy3 PRIVATE HYPHY3_COMBINED_DRIVER)
target_link_libraries(hyphy3 PRIVATE hyphy3_report)

# Standalone CLI Applications
add_executable(hyphy3-fel src/apps/hyphy_fel.cpp)
target_link_libraries(hyphy3-fel PRIVATE hyphy3_report)

add_executable(hyphy3-contrast-fel src/apps/hyphy_contrast_fel.cpp)
target_link_libraries(hyphy3-contrast-fel PRIVATE hyphy3_core)

add_executable(hyphy3-meme src/apps/hyphy_meme.cpp)
target_link_libraries(hyphy3-meme PRIVATE hyphy3_report)

add_executable(hyphy3-busted src/apps/hyphy_busted.cpp)
target_link_libraries(hyphy3-busted PRIVATE hyphy3_core)

add_executable(hyphy3-absrel src/apps/hyphy_absrel.cpp)
target_link_libraries(hyphy3-absrel PRIVATE hyphy3_core)

add_executable(hyphy3-relax src/apps/hyphy_relax.cpp)
target_link_libraries(hyphy3-relax PRIVATE hyphy3_core)

add_executable(hyphy3-fubar src/apps/hyphy_fubar.cpp)
target_link_libraries(hyphy3-fubar PRIVATE hyphy3_core)

add_executable(hyphy3-prime src/apps/hyphy_prime.cpp)
target_link_libraries(hyphy3-prime PRIVATE hyphy3_core)

add_executable(hyphy3-gard src/apps/hyphy_gard.cpp)
target_link_libraries(hyphy3-gard PRIVATE hyphy3_core)

add_executable(hyphy3-slac src/apps/hyphy_slac.cpp)
target_link_libraries(hyphy3-slac PRIVATE hyphy3_core)

add_executable(hyphy3-dfel src/apps/hyphy_dfel.cpp)
target_link_libraries(hyphy3-dfel PRIVATE hyphy3_report)

add_executable(hyphy3-ohmic src/apps/hyphy_ohmic.cpp)
target_link_libraries(hyphy3-ohmic PRIVATE hyphy3_report)

# Backward-compatible aliases
add_executable(hyphy-fel ALIAS hyphy3-fel)
add_executable(hyphy-contrast-fel ALIAS hyphy3-contrast-fel)
add_executable(hyphy-meme ALIAS hyphy3-meme)
add_executable(hyphy-busted ALIAS hyphy3-busted)
add_executable(hyphy-absrel ALIAS hyphy3-absrel)
add_executable(hyphy-relax ALIAS hyphy3-relax)
add_executable(hyphy-fubar ALIAS hyphy3-fubar)
add_executable(hyphy-prime ALIAS hyphy3-prime)
add_executable(hyphy-gard ALIAS hyphy3-gard)
add_executable(hyphy-slac ALIAS hyphy3-slac)
add_executable(hyphy-dfel ALIAS hyphy3-dfel)
add_executable(hyphy-ohmic ALIAS hyphy3-ohmic)

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

