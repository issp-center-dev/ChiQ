include(EnableGtests) #defined in ./cmake
set(test_src
  blockmatrixTest
  # bseTest
)

foreach(test ${test_src})
  add_gtest(${test})
endforeach(test)
