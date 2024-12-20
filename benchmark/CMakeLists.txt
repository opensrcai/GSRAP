
add_executable(bench_pnp_solver bench_pnp_solver.cc)
list(APPEND GSRAP_BENCHES bench_pnp_solver)

add_executable(bench_p3p_solver bench_p3p_solver.cc)
list(APPEND GSRAP_BENCHES bench_p3p_solver)

foreach(GSRAP_BENCH ${GSRAP_BENCHES})
  target_include_directories(${GSRAP_BENCH} PRIVATE ${PROJECT_SOURCE_DIR}/example)
  target_include_directories(${GSRAP_BENCH} PRIVATE ${PROJECT_SOURCE_DIR}/third_party/cereal/include)
  target_link_libraries(${GSRAP_BENCH} PRIVATE ${PROJECT_NAME}  cxxopts::cxxopts)
  target_compile_features(${GSRAP_BENCH} PRIVATE cxx_std_17)

  target_compile_options(${GSRAP_BENCH} PRIVATE 
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-global-constructors>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wall -Wextra -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-global-constructors>
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wredundant-decls -Wcast-align -Wmissing-declarations -Wmissing-include-dirs -Wswitch-enum -Wswitch-default -Winvalid-pch -Wredundant-decls -Wformat=2 -Wmissing-format-attribute -Wformat-nonliteral>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
  )

  if (${${PROJECT_PREFIX}_WARNINGS_AS_ERRORS})
    target_compile_options(${GSRAP_BENCH} PRIVATE 
      $<$<CXX_COMPILER_ID:Clang>:-Werror>
      $<$<CXX_COMPILER_ID:AppleClang>:-Werror>
      $<$<CXX_COMPILER_ID:GNU>:-Werror>
      $<$<CXX_COMPILER_ID:MSVC>:/wx>
    )
  endif (${${PROJECT_PREFIX}_WARNINGS_AS_ERRORS})
endforeach(GSRAP_BENCH ${GSRAP_BENCH})
