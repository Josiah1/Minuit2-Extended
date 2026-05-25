# Minuit2-Extended examples
# Built only when BUILD_EXAMPLES=ON

find_package(ROOT 6.0 REQUIRED COMPONENTS Core Hist MathCore Minuit Minuit2)
include(${ROOT_USE_FILE})

set(EXAMPLES
    basic_fit
    algorithms
    contour_scan
    profile_likelihood
    grid_scan
    error_budget
    multi_start
    state_snapshot
    iteration_trace
)

foreach(example ${EXAMPLES})
    add_executable(${example} ${example}.cc)
    target_link_libraries(${example} PRIVATE minuit2_extended ${ROOT_LIBRARIES})
endforeach()
