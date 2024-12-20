CC := clang
CXX := clang++

CMAKE_BUILD_RELEASE_DIR := cmake-build-release
CMAKE_BUILD_DEBUG_DIR := cmake-build-debug

CMAKE_COMMON_OPTIONS := -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DGSRAP_FIND_EIGEN=ON -DGSRAP_FIND_SPDLOG=ON -DGSRAP_FIND_THREADS=ON -DGSRAP_BUILD_TESTS=ON -DGSRAP_BUILD_EXAMPLES=ON -DGSRAP_BUILD_BENCHMARKS=ON -DGSRAP_ENABLE_NANOBIND=ON -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-linux-release -DVCPKG_MANIFEST_MODE=ON -DCMAKE_MODULE_PATH=$(shell dirname $(shell which python))/../lib/python$(shell python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')/site-packages/nanobind/cmake/ -DCMAKE_POSITION_INDEPENDENT_CODE=ON

PHONY_BASE := all release debug common cmake_configure_release cmake_configure_debug cmake_build_release cmake_build_debug weak_clean clean
.PHONY: $(PHONY_BASE)

release: cmake_configure_release cmake_build_release common
debug: cmake_configure_debug cmake_build_debug common

all: release debug

common:

# CMake Release Configure
cmake_configure_release: $(CMAKE_BUILD_RELEASE_DIR)/Makefile
$(CMAKE_BUILD_RELEASE_DIR)/Makefile:
	CC=$(CC) CXX=$(CXX) cmake -H. -B$(CMAKE_BUILD_RELEASE_DIR) -DCMAKE_BUILD_TYPE=Release $(CMAKE_COMMON_OPTIONS)

# CMake Debug Configure
cmake_configure_debug: $(CMAKE_BUILD_DEBUG_DIR)/Makefile
$(CMAKE_BUILD_DEBUG_DIR)/Makefile:
	CC=$(CC) CXX=$(CXX) cmake -H. -B$(CMAKE_BUILD_DEBUG_DIR) -DCMAKE_BUILD_TYPE=Debug $(CMAKE_COMMON_OPTIONS)

# CMake Release Build
cmake_build_release: cmake_configure_release
	cd $(CMAKE_BUILD_RELEASE_DIR) && $(MAKE)

# CMake Debug Build
cmake_build_debug: cmake_configure_debug
	cd $(CMAKE_BUILD_DEBUG_DIR) && $(MAKE)

weak_clean:
	if [ -f $(CMAKE_BUILD_RELEASE_DIR)/Makefile ]; then cd $(CMAKE_BUILD_RELEASE_DIR) && $(MAKE) clean && cd .. ; fi
	if [ -f $(CMAKE_BUILD_DEBUG_DIR)/Makefile ]; then cd $(CMAKE_BUILD_DEBUG_DIR) && $(MAKE) clean && cd .. ; fi

	# if [ -f $(CMAKE_BUILD_RELEASE_DIR)/Makefile ]; then cmake --build $(CMAKE_BUILD_RELEASE_DIR) --target clean ; fi
	# if [ -f $(CMAKE_BUILD_DEBUG_DIR)/Makefile ]; then cmake --build $(CMAKE_BUILD_DEBUG_DIR) --target clean ; fi

clean: weak_clean
	@rm -rf $(CMAKE_BUILD_RELEASE_DIR) $(CMAKE_BUILD_DEBUG_DIR)
