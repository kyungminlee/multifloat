project(multifloats LANGUAGES CXX)

include(GNUInstallDirs)

add_library(multifloats STATIC multifloats_math.cc)
target_include_directories(multifloats
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# Hide all C++ symbols; only extern "C" dd_* functions are exported.
set_target_properties(multifloats PROPERTIES CXX_VISIBILITY_PRESET hidden)
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(multifloats PRIVATE -O3 -flto -ffat-lto-objects)
    else()
        target_compile_options(multifloats PRIVATE -O3 -flto)
    endif()
else()
    target_compile_options(multifloats PRIVATE -O2)
endif()

# Install: library + header + export targets + Config.cmake
install(TARGETS multifloats
    EXPORT multifloatsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES multifloats.hh
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT multifloatsTargets
    FILE multifloatsTargets.cmake
    NAMESPACE multifloats::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloats)

include(CMakePackageConfigHelpers)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfig.cmake"
"include(\"\${CMAKE_CURRENT_LIST_DIR}/multifloatsTargets.cmake\")
")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloats)
