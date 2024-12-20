#!/bin/bash

set -e

find . -type f -regextype posix-basic -regex "./\(src\|example\|test\|benchmark\)/.*\.cc" | xargs -L 1 -P 24 clang-tidy --extra-arg=-std=c++17 --fix --fix-errors
