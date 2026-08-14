set(BUILD_SHARED_LIBS OFF)

add_subdirectory("toml-f" EXCLUDE_FROM_ALL)
add_subdirectory("mt_stream_f90-1.11")

list(
  APPEND lib-deps
    "toml-f-lib"
    "mt_stream-lib"
)

set(lib-deps "${lib-deps}" PARENT_SCOPE)
