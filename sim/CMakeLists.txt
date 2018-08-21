# The main process executable, for debugging and playing with the code.
add_executable(template template.cpp)
target_link_libraries(template ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})
