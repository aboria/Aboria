file(GLOB Tyche_SOURCES "*.cpp")
file(GLOB Tyche_HEADERS "*.h")
include_directories(${Tyche_INCLUDE_DIRECTORIES})
message(STATUS ${Tyche_INCLUDE_DIRECTORIES})

add_library(Tyche SHARED ${Tyche_SOURCES} "../smoldyn/rxnparam.c")
TARGET_LINK_LIBRARIES(Tyche ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${Boost_TIMER_LIBRARY})

install (TARGETS Tyche DESTINATION lib)
install (FILES ${Tyche_HEADERS} DESTINATION include)
install (FILES ../smoldyn/rxnparam.h DESTINATION include)
install (DIRECTORY ${Tyche_SOURCE_DIR}/eigen3/ DESTINATION include)
