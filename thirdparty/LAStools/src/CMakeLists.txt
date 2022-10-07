if (!MSVC)
add_compile_options(-Wno-deprecated -Wno-write-strings -Wno-unused-result)
endif()
add_definitions(-DNDEBUG )

include_directories(../LASzip/src)
include_directories(../LASlib/inc)

set(GEOPROJECTION_TARGETS
	las2las
)

set(ALL_TARGETS ${GEOPROJECTION_TARGETS})

foreach(TARGET ${GEOPROJECTION_TARGETS})
	add_executable(${TARGET} ${TARGET}.cpp geoprojectionconverter.cpp) 
endforeach(TARGET)

foreach(TARGET ${ALL_TARGETS})
	target_link_libraries(${TARGET} LASlib)
	set_target_properties(${TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools/LAStools)
	set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME city4cfd_${TARGET})
	install(TARGETS ${TARGET} RUNTIME DESTINATION bin)
endforeach(TARGET)
