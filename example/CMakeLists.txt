set(EXECUTABLE_TARGETS "")

add_executable(example_essential_solver example_essential_solver.cc)
list(APPEND EXECUTABLE_TARGETS example_essential_solver)

add_executable(example_pnp_solver example_pnp_solver.cc)
list(APPEND EXECUTABLE_TARGETS example_pnp_solver)

add_executable(example_p3p_solver example_p3p_solver.cc)
list(APPEND EXECUTABLE_TARGETS example_p3p_solver)

add_executable(example_sim3_solver example_sim3_solver.cc)
list(APPEND EXECUTABLE_TARGETS example_sim3_solver)


foreach(EXECUTABLE_TARGET IN LISTS EXECUTABLE_TARGETS)
    set_target_properties(${EXECUTABLE_TARGET} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

    target_compile_features(${EXECUTABLE_TARGET} PRIVATE cxx_std_17)

    target_compile_options(${EXECUTABLE_TARGET} PRIVATE 
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-global-constructors>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wall -Wextra -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-global-constructors>
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wredundant-decls -Wcast-align -Wmissing-declarations -Wmissing-include-dirs -Wswitch-enum -Wswitch-default -Winvalid-pch -Wredundant-decls -Wformat=2 -Wmissing-format-attribute -Wformat-nonliteral>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    )    

    if (${${PROJECT_PREFIX}_WARNINGS_AS_ERRORS})
      target_compile_options(${EXECUTABLE_TARGET} PRIVATE 
        $<$<CXX_COMPILER_ID:Clang>:-Werror>
        $<$<CXX_COMPILER_ID:AppleClang>:-Werror>
        $<$<CXX_COMPILER_ID:GNU>:-Werror>
        $<$<CXX_COMPILER_ID:MSVC>:/wx>
      )
    endif (${${PROJECT_PREFIX}_WARNINGS_AS_ERRORS})

    target_include_directories(${EXECUTABLE_TARGET} PRIVATE ${PROJECT_SOURCE_DIR}/inc)
    target_include_directories(${EXECUTABLE_TARGET} PRIVATE ${PROJECT_SOURCE_DIR}/example)

    target_link_libraries(${EXECUTABLE_TARGET} PRIVATE ${PROJECT_NAME})

endforeach()
