#  Interpolate a field that is cubic in R and Z from an mgrid table with becoil.
add_executable (xbecoil_test ${CMAKE_CURRENT_SOURCE_DIR}/test_becoil.f90)
target_link_libraries (xbecoil_test PUBLIC vmec)

add_test (NAME    vmec_becoil_cubic_test
          COMMAND $<TARGET_FILE:xbecoil_test>)
