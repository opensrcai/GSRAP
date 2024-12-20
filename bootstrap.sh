#!/bin/bash

rm -rf cmake-build-debug
rm -rf cmake-build-release

CC=clang CXX=clang++ cmake -Bcmake-build-debug -H. \
                           -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
                           -DCMAKE_BUILD_TYPE=Debug \
                           -DGSRAP_FIND_EIGEN=ON \
                           -DGSRAP_FIND_SPDLOG=ON \
                           -DGSRAP_FIND_THREADS=ON \
                           -DGSRAP_BUILD_TESTS=ON \
                           -DGSRAP_BUILD_EXAMPLES=ON \
                           -DGSRAP_BUILD_BENCHMARKS=ON \
                           -DGSRAP_BUILD_SHARE=ON \
                           -DGSRAP_OPT_FOR_CPU=OFF \
                           -DGSRAP_SHOW_OPT_RESULT=OFF \
                           -DGSRAP_WARNINGS_AS_ERRORS=ON \
                           -DGSRAP_USE_SANITIZER=ON

mv cmake-build-debug/compile_commands.json .

CC=clang CXX=clang++ cmake -Bcmake-build-release -H. \
                           -DCMAKE_EXPORT_COMPILE_COMMANDS=0 \
                           -DCMAKE_BUILD_TYPE=Release \
                           -DGSRAP_FIND_EIGEN=ON \
                           -DGSRAP_FIND_SPDLOG=ON \
                           -DGSRAP_FIND_THREADS=ON \
                           -DGSRAP_BUILD_TESTS=ON \
                           -DGSRAP_BUILD_EXAMPLES=ON \
                           -DGSRAP_BUILD_BENCHMARKS=ON \
                           -DGSRAP_BUILD_SHARE=ON \
                           -DGSRAP_OPT_FOR_CPU=ON \
                           -DGSRAP_SHOW_OPT_RESULT=OFF \
                           -DGSRAP_WARNINGS_AS_ERRORS=ON \
                           -DGSRAP_USE_SANITIZER=OFF
