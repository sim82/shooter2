project(shooter2)

#add_executable(shooter2 main.cpp)
#add_executable(rad rad.cpp)
#add_executable(clangl clangl.cpp)
add_executable( ortho ortho.cpp player_bits.cpp scene_bits.cpp crystal_bits.cpp rad_core.cpp gl_bits.cpp )


SET(ClanLib_DIR  ${CMAKE_CURRENT_SOURCE_DIR})
SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall -std=c++0x -march=native" )
# SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall -std=c++0x -march=core2" )

SET( ClanLib_MAJOR_VERSION 2 )
SET( ClanLib_MINOR_VERSION 3 )

FIND_PACKAGE(ClanLib REQUIRED)

INCLUDE_DIRECTORIES(${ClanLib_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(ortho ${ClanLib_Core_LIBRARY} ${ClanLib_App_LIBRARY} ${ClanLib_Sound_LIBRARY} ${ClanLib_Display_LIBRARY} ${ClanLib_GL_LIBRARY}  ${ClanLib_GUI_LIBRARY}  -lclan23CSSLayout -lGL -lX11 -lfreetype -lfontconfig -pthread OpenCL)
#TARGET_LINK_LIBRARIES(clangl ${ClanLib_Core_LIBRARY} ${ClanLib_App_LIBRARY} ${ClanLib_Sound_LIBRARY} ${ClanLib_Display_LIBRARY} ${ClanLib_GL_LIBRARY}  ${ClanLib_GUI_LIBRARY}  -lclan23CSSLayout -lGL -lX11 -lfreetype -lfontconfig -pthread)
#TARGET_LINK_LIBRARIES(shooter2 ${ClanLib_Core_LIBRARY} ${ClanLib_App_LIBRARY} ${ClanLib_Sound_LIBRARY} ${ClanLib_Display_LIBRARY} ${ClanLib_GL_LIBRARY} -lX11 -lfreetype -lfontconfig -pthread)
#TARGET_LINK_LIBRARIES(rad ${ClanLib_Core_LIBRARY} ${ClanLib_App_LIBRARY} ${ClanLib_Sound_LIBRARY} ${ClanLib_Display_LIBRARY} ${ClanLib_GL_LIBRARY}  -lX11 -lfreetype -lfontconfig -pthread)
