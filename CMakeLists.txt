# Minuit2-Extended: template-based Minuit2/Minuit1 wrapper library
# Header-only INTERFACE library — no separate compilation needed.

cmake_minimum_required(VERSION 3.15)

project(Minuit2Extended
    VERSION 2.0.0
    DESCRIPTION "Template-based Minuit2/Minuit1 minimization wrapper"
    LANGUAGES CXX
)

# When used as a subdirectory, ROOT is found by the parent project.
# When built standalone (examples), find ROOT here.
if(NOT ROOT_FOUND)
    find_package(ROOT 6.0 COMPONENTS Core Hist MathCore Minuit Minuit2)
endif()

add_library(minuit2_extended INTERFACE)

target_include_directories(minuit2_extended INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(minuit2_extended INTERFACE cxx_std_17)

# Propagate ROOT include dirs and libraries when available
if(ROOT_FOUND)
    target_include_directories(minuit2_extended INTERFACE ${ROOT_INCLUDE_DIRS})
    target_link_libraries(minuit2_extended INTERFACE
        ROOT::Minuit2
        ROOT::MathCore
        ROOT::Hist
        ROOT::Core
    )
endif()

# Install
install(TARGETS minuit2_extended EXPORT Minuit2ExtendedTargets)
install(DIRECTORY include/ DESTINATION include)
install(EXPORT Minuit2ExtendedTargets
    FILE Minuit2ExtendedConfig.cmake
    NAMESPACE Minuit2Extended::
    DESTINATION lib/cmake/Minuit2Extended
)

# Examples (opt-in)
option(BUILD_EXAMPLES "Build Minuit2-Extended examples" OFF)
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
