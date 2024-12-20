set(${PROJECT_PREFIX}_SOURCE 
  ${CMAKE_CURRENT_SOURCE_DIR}/p3p_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/p3p_solver/compute_poses.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/pnp_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/essential_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sim3_solver.cc
  )

set(${PROJECT_PREFIX}_THIRD_PARTY_SOURCE 
  ${PROJECT_SOURCE_DIR}/third_party/opengv/src/absolute_pose/Epnp.cpp
  ${PROJECT_SOURCE_DIR}/third_party/opengv/src/math/Sturm.cpp
  ${PROJECT_SOURCE_DIR}/third_party/opengv/src/relative_pose/modules/main.cpp
  ${PROJECT_SOURCE_DIR}/third_party/opengv/src/relative_pose/modules/fivept_nister/modules.cpp)

if (${${PROJECT_PREFIX}_BUILD_SHARE})
  add_library(${PROJECT_NAME} SHARED
    ${${PROJECT_PREFIX}_SOURCE} ${${PROJECT_PREFIX}_THIRD_PARTY_SOURCE})
else (${${PROJECT_PREFIX}_BUILD_SHARE})
  add_library(${PROJECT_NAME} STATIC
    ${${PROJECT_PREFIX}_SOURCE} ${${PROJECT_PREFIX}_THIRD_PARTY_SOURCE})
endif (${${PROJECT_PREFIX}_BUILD_SHARE})

set_source_files_properties(${${PROJECT_PREFIX}_SOURCE} PROPERTIES
  COMPILE_FLAGS "$<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded>$<$<CXX_COMPILER_ID:AppleClang>:-Wall -Wextra -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded>$<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wredundant-decls -Wcast-align -Wmissing-declarations -Wmissing-include-dirs -Wswitch-enum -Wswitch-default -Winvalid-pch -Wredundant-decls -Wformat=2 -Wmissing-format-attribute -Wformat-nonliteral>$<$<CXX_COMPILER_ID:MSVC>:/W4>"
  )

if (${${PROJECT_PREFIX}_WARNINGS_AS_ERRORS})
  set_source_files_properties(${${PROJECT_PREFIX}_SOURCE} PROPERTIES
    COMPILE_FLAGS "$<$<CXX_COMPILER_ID:Clang>:-Werror>$<$<CXX_COMPILER_ID:AppleClang>:-Werror>$<$<CXX_COMPILER_ID:GNU>:-Werror>$<$<CXX_COMPILER_ID:MSVC>:/wx>"
    )
endif (${${PROJECT_PREFIX}_WARNINGS_AS_ERRORS})

set_source_files_properties(${${PROJECT_PREFIX}_THIRD_PARTY_SOURCE} PROPERTIES
  COMPILE_FLAGS "$<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra>$<$<CXX_COMPILER_ID:AppleClang>:-Wall -Wextra>$<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra>$<$<CXX_COMPILER_ID:MSVC>:/W3>")

if (${${PROJECT_PREFIX}_OPT_FOR_CPU})
  target_compile_options(${PROJECT_NAME} PRIVATE $<$<CXX_COMPILER_ID:Clang>:-march=native>$<$<CXX_COMPILER_ID:AppleClang>:-march=native>$<$<CXX_COMPILER_ID:GNU>:-march=native>$<$<CXX_COMPILER_ID:MSVC>:/arch:AVX2>)
endif (${${PROJECT_PREFIX}_OPT_FOR_CPU})

message("")
if (${${PROJECT_PREFIX}_SHOW_OPT_RESULT})
  target_compile_options(${PROJECT_NAME} PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Rpass=loop-vectorize>$<$<CXX_COMPILER_ID:AppleClang>:-Rpass=loop-vectorize>$<$<CXX_COMPILER_ID:GNU>:-fopt-info-vec-optimized>$<$<CXX_COMPILER_ID:MSVC>:>)
endif (${${PROJECT_PREFIX}_SHOW_OPT_RESULT})

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/inc)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/third_party/opengv/include)

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} PRIVATE spdlog::spdlog Threads::Threads)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

# target_compile_options(${PROJECT_NAME} PUBLIC -DGSRAP_LESS_EIGEN_3_4)
