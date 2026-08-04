add_subdirectory(common)
add_subdirectory(reconstruction)

add_library(roofer-core)
target_link_libraries(roofer-core PUBLIC ptinpoly common reconstruction)

# The logger (and its fmt dependency) is only needed by apps, tests, and the
# doc-helper - nothing in common/reconstruction uses it.
if(RF_BUILD_APPS OR RF_BUILD_TESTING OR RF_BUILD_DOC_HELPER)
  add_subdirectory(logger)
  target_link_libraries(roofer-core PUBLIC logger)
endif()
