# End-to-end check that `hyphy3 <method> --html` and `hyphy3 report` produce standalone reports.
# Run as: cmake -DHYPHY3=<binary> -DDATA_DIR=<benchmarks/data> -DWORK_DIR=<dir> -P check_html_report.cmake

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

function(run)
    execute_process(COMMAND ${ARGN} WORKING_DIRECTORY "${WORK_DIR}"
                    RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "Command failed (${rc}): ${ARGN}\n${out}\n${err}")
    endif()
endfunction()

function(check_report html tag)
    if(NOT EXISTS "${html}")
        message(FATAL_ERROR "Report not written: ${html}")
    endif()
    file(READ "${html}" content)
    # The header takes the producer chip the viewer names on load
    foreach(needle "<${tag}" "window.HyPhyFELViewer" ".hyphy-fel-container" "const payload = {"
                   "<div class=\"app-nav-meta\">" "addEventListener('data-loaded', showFormat)")
        string(FIND "${content}" "${needle}" pos)
        if(pos EQUAL -1)
            message(FATAL_ERROR "${html} is missing '${needle}'")
        endif()
    endforeach()
endfunction()

# 1. Report written as part of an analysis run (small alignment keeps the test fast)
run("${HYPHY3}" fel --alignment "${DATA_DIR}/cd2_reduced.fna" --tree "${DATA_DIR}/cd2_reduced.nwk"
    --threads 2 --no-progress --output fel.json --html fel.html)
check_report("${WORK_DIR}/fel.html" hyphy-fel-viewer)

# 2. Report from an existing result file, including the legacy Datamonkey format
run("${HYPHY3}" fel --alignment "${DATA_DIR}/cd2_reduced.fna" --tree "${DATA_DIR}/cd2_reduced.nwk"
    --threads 2 --no-progress --legacy-json --output fel.legacy.json)
run("${HYPHY3}" report fel.legacy.json)
check_report("${WORK_DIR}/fel.legacy.html" hyphy-fel-viewer)

# 3. FUBAR reports mount the neutral viewer element, from a run and from its result file
run("${HYPHY3}" fubar --alignment "${DATA_DIR}/cd2_reduced.fna" --tree "${DATA_DIR}/cd2_reduced.nwk"
    --threads 2 --no-progress --output fubar.json --html fubar.html)
check_report("${WORK_DIR}/fubar.html" hyphy-viewer)
run("${HYPHY3}" report fubar.json --html fubar.report.html)
check_report("${WORK_DIR}/fubar.report.html" hyphy-viewer)

# 4. Analyses without a viewer are refused with a clear error rather than mis-rendered
file(WRITE "${WORK_DIR}/busted.json" "{\"analysis\": {\"id\": \"busted\"}}")
execute_process(COMMAND "${HYPHY3}" report busted.json --html busted.html WORKING_DIRECTORY "${WORK_DIR}"
                RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_QUIET)
if(rc EQUAL 0 OR EXISTS "${WORK_DIR}/busted.html")
    message(FATAL_ERROR "hyphy3 report should refuse results without a viewer")
endif()
string(FIND "${err}" "no HTML viewer for 'busted'" pos)
if(pos EQUAL -1)
    message(FATAL_ERROR "Unexpected refusal message: ${err}")
endif()

message(STATUS "HTML report checks passed")
