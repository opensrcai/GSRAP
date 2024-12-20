find_package(GTest CONFIG REQUIRED)
include(GoogleTest)

add_executable(test_essential_solver test_essential_solver.cc)
list(APPEND GSRAP_TESTS test_essential_solver)

add_executable(test_pnp_solver test_pnp_solver.cc)
list(APPEND GSRAP_TESTS test_pnp_solver)

add_executable(test_p3p_solver test_p3p_solver.cc)
list(APPEND GSRAP_TESTS test_p3p_solver)

add_executable(test_sim3_solver test_sim3_solver.cc)
list(APPEND GSRAP_TESTS test_sim3_solver)

foreach(GSRAP_TEST ${GSRAP_TESTS})
  target_include_directories(${GSRAP_TEST} PRIVATE ${PROJECT_SOURCE_DIR}/example)
  target_link_libraries(${GSRAP_TEST} PRIVATE ${PROJECT_NAME} GTest::gmock GTest::gtest GTest::gmock_main GTest::gtest_main)
  target_compile_features(${GSRAP_TEST} PRIVATE cxx_std_17)

  target_compile_options(${GSRAP_TEST} PRIVATE 
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-global-constructors>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wall -Wextra -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-global-constructors>
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wredundant-decls -Wcast-align -Wmissing-declarations -Wmissing-include-dirs -Wswitch-enum -Wswitch-default -Winvalid-pch -Wredundant-decls -Wformat=2 -Wmissing-format-attribute -Wformat-nonliteral>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
  )

  if (${${PROJECT_PREFIX}_WARNINGS_AS_ERRORS})
    target_compile_options(${GSRAP_TEST} PRIVATE 
      $<$<CXX_COMPILER_ID:Clang>:-Werror>
      $<$<CXX_COMPILER_ID:AppleClang>:-Werror>
      $<$<CXX_COMPILER_ID:GNU>:-Werror>
      $<$<CXX_COMPILER_ID:MSVC>:/wx>
    )
  endif (${${PROJECT_PREFIX}_WARNINGS_AS_ERRORS})

  gtest_discover_tests(${GSRAP_TEST}) # https://mickey-happygolucky.hatenablog.com/entry/2021/02/04/005424
endforeach(GSRAP_TEST ${GSRAP_TESTS})
