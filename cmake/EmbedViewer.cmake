# Bundles the HTML results viewer (viz/fel-viewer) into a C++ header so hyphy3 can write
# self-contained HTML reports without Node, Python or any files next to the binary.
#
# Run in script mode:
#   cmake -DVIEWER_DIR=<viz/fel-viewer> -DOUTPUT=<header> [-DSTAMP=<file>] -P EmbedViewer.cmake
#
# The viewer is plain ES modules with only static, relative, named imports. Bundling is
# therefore: resolve imports depth-first from src/index.js, emit each module once in
# dependency order, drop the import lines and `export` keywords, and wrap everything in
# one function scope that publishes the entry module's exports as window.HyPhyFELViewer
# (the same global the former UMD build exposed) and under the neutral window.HyPhyViewer.

cmake_minimum_required(VERSION 3.20)

if(NOT VIEWER_DIR OR NOT OUTPUT)
    message(FATAL_ERROR "EmbedViewer.cmake needs -DVIEWER_DIR=... and -DOUTPUT=...")
endif()

set(_ordered "")
set(_visiting "")

function(_visit module_path)
    get_filename_component(module_path "${module_path}" REALPATH)
    if(module_path IN_LIST _ordered)
        return()
    endif()
    if(module_path IN_LIST _visiting)
        message(FATAL_ERROR "Circular import involving ${module_path}")
    endif()
    list(APPEND _visiting "${module_path}")
    set(_visiting "${_visiting}" PARENT_SCOPE)

    if(NOT EXISTS "${module_path}")
        message(FATAL_ERROR "Viewer module not found: ${module_path}")
    endif()
    file(STRINGS "${module_path}" _imports REGEX "^import ")
    get_filename_component(_dir "${module_path}" DIRECTORY)
    foreach(_line IN LISTS _imports)
        if(NOT _line MATCHES "from '(\\./|\\.\\./)[^']+\\.js'")
            message(FATAL_ERROR "Unsupported import in ${module_path}: ${_line}\n"
                                "Only static named imports of relative .js modules can be embedded.")
        endif()
        string(REGEX REPLACE "^.*from '([^']+)'.*$" "\\1" _rel "${_line}")
        _visit("${_dir}/${_rel}")
        set(_ordered "${_ordered}" PARENT_SCOPE)
    endforeach()

    list(APPEND _ordered "${module_path}")
    set(_ordered "${_ordered}" PARENT_SCOPE)
endfunction()

set(_entry "${VIEWER_DIR}/src/index.js")
_visit("${_entry}")

# Names exported by the entry module become the public API
file(STRINGS "${_entry}" _entry_exports REGEX "^export (class|function|const|let) ")
set(_api "")
foreach(_line IN LISTS _entry_exports)
    string(REGEX REPLACE "^export (class|function|const|let) ([A-Za-z_$][A-Za-z0-9_$]*).*$" "\\2" _name "${_line}")
    list(APPEND _api "${_name}")
endforeach()
list(JOIN _api ", " _api_list)

# All modules share one scope after bundling, so top-level names must be unique across modules
set(_declared "")
foreach(_module IN LISTS _ordered)
    file(STRINGS "${_module}" _decls REGEX "^(export )?(async )?(class|function|const|let|var) [A-Za-z_$][A-Za-z0-9_$]*")
    foreach(_decl IN LISTS _decls)
        string(REGEX REPLACE "^(export )?(async )?(class|function|const|let|var) ([A-Za-z_$][A-Za-z0-9_$]*).*$" "\\4" _name "${_decl}")
        if(_name IN_LIST _declared)
            message(FATAL_ERROR "Top-level name '${_name}' is declared in more than one viewer module "
                                "(again in ${_module}). Bundled modules share one scope; rename one of them.")
        endif()
        list(APPEND _declared "${_name}")
    endforeach()
endforeach()

set(_js "(function () {\n'use strict';\n")
foreach(_module IN LISTS _ordered)
    file(READ "${_module}" _src)
    file(RELATIVE_PATH _rel "${VIEWER_DIR}" "${_module}")
    # Drop import statements (one per line in this codebase) and strip `export` keywords
    string(REGEX REPLACE "(^|\n)import [^\n]*" "\\1" _src "${_src}")
    string(REGEX REPLACE "(^|\n)export (class|function|const|let|async function) " "\\1\\2 " _src "${_src}")
    if(_src MATCHES "(^|\n)export ")
        message(FATAL_ERROR "Unsupported export form in ${_module}; use `export class|function|const|let`.")
    endif()
    string(APPEND _js "\n// ---- ${_rel} ----\n${_src}\n")
endforeach()
string(APPEND _js "\nwindow.HyPhyFELViewer = window.HyPhyViewer = { ${_api_list} };\n})();\n")

# Stylesheets are linked from the dev page rather than imported, so embed them all
file(GLOB _css_files "${VIEWER_DIR}/src/styles/*.css")
list(SORT _css_files)
set(_css "")
foreach(_file IN LISTS _css_files)
    file(READ "${_file}" _content)
    string(APPEND _css "${_content}\n")
endforeach()

# Emit as chunks of raw string literals: some compilers cap a single literal at 64 KB
set(_delim "HYV3")
foreach(_text IN ITEMS _js _css)
    string(FIND "${${_text}}" ")${_delim}\"" _clash)
    if(NOT _clash EQUAL -1)
        message(FATAL_ERROR "Viewer source contains the raw-string delimiter )${_delim}\"")
    endif()
endforeach()

function(_chunks text out_var)
    string(LENGTH "${text}" _len)
    set(_pos 0)
    set(_result "")
    while(_pos LESS _len)
        string(SUBSTRING "${text}" ${_pos} 16000 _piece)
        string(APPEND _result "    R\"${_delim}(${_piece})${_delim}\",\n")
        math(EXPR _pos "${_pos} + 16000")
    endwhile()
    set(${out_var} "${_result}" PARENT_SCOPE)
endfunction()

_chunks("${_js}" _js_chunks)
_chunks("${_css}" _css_chunks)
list(LENGTH _ordered _n_modules)

set(_header "// Generated by cmake/EmbedViewer.cmake from viz/fel-viewer. Do not edit.
#pragma once

#include <string>

namespace hyphy::report::assets {

// ${_n_modules} modules from viz/fel-viewer/src, bundled in dependency order
inline const std::string& viewer_js() {
    static const char* const chunks[] = {
${_js_chunks}    };
    static const std::string joined = [] {
        std::string s;
        for (const char* c : chunks) s += c;
        return s;
    }();
    return joined;
}

inline const std::string& viewer_css() {
    static const char* const chunks[] = {
${_css_chunks}    \"\"
    };
    static const std::string joined = [] {
        std::string s;
        for (const char* c : chunks) s += c;
        return s;
    }();
    return joined;
}

}  // namespace hyphy::report::assets
")

# Rewrite the header only when its content changes, so dependent targets don't rebuild needlessly;
# the stamp file (the build rule's real output) is always refreshed.
set(_changed TRUE)
if(EXISTS "${OUTPUT}")
    file(READ "${OUTPUT}" _old)
    if(_old STREQUAL _header)
        set(_changed FALSE)
    endif()
endif()
if(_changed)
    file(WRITE "${OUTPUT}" "${_header}")
    message(STATUS "Embedded ${_n_modules} viewer modules into ${OUTPUT}")
endif()
if(STAMP)
    file(TOUCH "${STAMP}")
endif()
