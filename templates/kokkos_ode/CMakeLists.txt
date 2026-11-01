# ABOUTME: CMake build file for JAFF-generated chemistry solver
# ABOUTME: Links against header-only integrators (VODE)

cmake_minimum_required(VERSION 3.16)
project(ChemistryODE LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Integration with header-only integrators library
# Default behavior: fetch from GitHub via FetchContent
option(USE_SYSTEM_INTEGRATORS "Use an installed integrators package (find_package)" OFF)
set(INTEGRATORS_GIT_REPOSITORY "https://github.com/BenWibking/header-only-integrators.git" CACHE STRING "Git repository for header-only-integrators")
set(INTEGRATORS_GIT_TAG "development" CACHE STRING "Git tag/branch/commit for header-only-integrators")

if(USE_SYSTEM_INTEGRATORS)
    find_package(integrators REQUIRED)
else()
    include(FetchContent)
    # Ensure dependency does not build its tests/examples
    set(_HOI_HAD_BUILD_TESTS FALSE)
    if(DEFINED BUILD_TESTS)
        set(_HOI_HAD_BUILD_TESTS TRUE)
        set(_HOI_SAVED_BUILD_TESTS "${BUILD_TESTS}")
    endif()
    set(_HOI_HAD_BUILD_EXAMPLES FALSE)
    if(DEFINED BUILD_EXAMPLES)
        set(_HOI_HAD_BUILD_EXAMPLES TRUE)
        set(_HOI_SAVED_BUILD_EXAMPLES "${BUILD_EXAMPLES}")
    endif()
    set(BUILD_TESTS OFF CACHE BOOL "Disable tests in header-only-integrators" FORCE)
    set(BUILD_EXAMPLES OFF CACHE BOOL "Disable examples in header-only-integrators" FORCE)
    FetchContent_Declare(header_only_integrators
        GIT_REPOSITORY ${INTEGRATORS_GIT_REPOSITORY}
        GIT_TAG        ${INTEGRATORS_GIT_TAG}
    )
    # This will add the 'integrators' target provided by the fetched project
    FetchContent_MakeAvailable(header_only_integrators)

    # Restore caller cache variables to avoid affecting the parent project
    if(_HOI_HAD_BUILD_TESTS)
        set(BUILD_TESTS "${_HOI_SAVED_BUILD_TESTS}" CACHE BOOL "" FORCE)
    else()
        unset(BUILD_TESTS CACHE)
    endif()
    if(_HOI_HAD_BUILD_EXAMPLES)
        set(BUILD_EXAMPLES "${_HOI_SAVED_BUILD_EXAMPLES}" CACHE BOOL "" FORCE)
    else()
        unset(BUILD_EXAMPLES CACHE)
    endif()
endif()

# Add main executable
add_executable(chemistry_ode chemistry_ode.cpp)

# Link libraries
target_link_libraries(chemistry_ode PRIVATE integrators)
