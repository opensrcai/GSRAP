/*
MIT License

Copyright (c) 2021 Okawa Kai

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#ifdef __clang__

// clang-format off
#define GSRAP_IGNORE_STRICT_WARNING_PUSH \
  _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Weverything\"")


#elif defined(__GNUC__) || defined(__GNUG__)

#define GSRAP_IGNORE_STRICT_WARNING_PUSH \
  _Pragma("GCC diagnostic ignored \"-Wshadow\"") \
  _Pragma("GCC diagnostic ignored \"-Wredundant-decls\"") \
  _Pragma("GCC diagnostic ignored \"-Wcast-align\"") \
  _Pragma("GCC diagnostic ignored \"-Wmissing-declarations\"") \
  _Pragma("GCC diagnostic ignored \"-Wmissing-include-dirs\"") \
  _Pragma("GCC diagnostic ignored \"-Wswitch-enum\"") \
  _Pragma("GCC diagnostic ignored \"-Wswitch-default\"") \
  _Pragma("GCC diagnostic ignored \"-Winvalid-pch\"") \
  _Pragma("GCC diagnostic ignored \"-Wredundant-decls\"") \
  _Pragma("GCC diagnostic ignored \"-Wformat=2\"") \
  _Pragma("GCC diagnostic ignored \"-Wmissing-format-attribute\"") \
  _Pragma("GCC diagnostic ignored \"-Wformat-nonliteral\"")

#elif defined(_MSC_VER)

#define GSRAP_IGNORE_STRICT_WARNING_PUSH \
  __pragma(warning(push)) \
  __pragma(warning(disable : 4706))

#else

#define GSRAP_IGNORE_STRICT_WARNING_PUSH

#endif
// clang-format on

#ifdef __clang__

#define GSRAP_IGNORE_STRICT_WARNING_POP _Pragma("clang diagnostic pop")

#elif defined(__GNUC__) || defined(__GNUG__)

#define GSRAP_IGNORE_STRICT_WARNING_POP _Pragma("GCC diagnostic pop")

#elif defined(_MSC_VER)

#define GSRAP_IGNORE_STRICT_WARNING_POP __pragma(warning(pop))

#else

#define GSRAP_IGNORE_STRICT_WARNING_POP

#endif
