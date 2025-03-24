cmake_minimum_required(VERSION 2.8.12)
include(EnableGtests) #defined in ./cmake
set(test_src
  blockmatrixTest
  # bseTest
)

foreach(test ${test_src})
  add_gtest(${test})
endforeach(test)
