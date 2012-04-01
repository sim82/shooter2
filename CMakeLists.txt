project(shooter2)

#add_executable(shooter2 main.cpp)
#add_executable(rad rad.cpp)
#add_executable(clangl clangl.cpp)
add_executable( ortho ortho.cpp player_bits.cpp scene_bits.cpp crystal_bits.cpp rad_core.cpp gl_bits.cpp )

IF(WIN32)
set( CLANLIB_ROOT E:\\src\\clanlib )
include_directories( ${include_directories} ${CLANLIB_ROOT}\\include)
LINK_DIRECTORIES( ${LINK_DIRECTORIES} ${CLANLIB_ROOT}\\lib\\x86 )
set(ClanLib_Core_LIBRARY clanCore-static-mt)
set(ClanLib_Core_LIBRARY clanApp-static-mt)
set(ClanLib_Core_LIBRARY clanSound-static-mt)
set(ClanLib_Core_LIBRARY clanDisplay-static-mt)
set(ClanLib_Core_LIBRARY clanGL-static-mt)
set(ClanLib_Core_LIBRARY clanGUI-static-mt)

include_directories( ${include_directories} E:\\src\\boost_1_48_0)
LINK_DIRECTORIES( ${LINK_DIRECTORIES} E:\\src\\boost_1_48_0\\stage\\lib)

set( BOOST_LIBS )
set( SYSDEP_LIBS )
ELSE()
SET( ClanLib_MAJOR_VERSION 2 )
SET( ClanLib_MINOR_VERSION 3 )

FIND_PACKAGE(ClanLib REQUIRED)

ENDIF()

SET(ClanLib_DIR  ${CMAKE_CURRENT_SOURCE_DIR})
SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall -std=c++0x" )


INCLUDE_DIRECTORIES(${ClanLib_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(ortho ${ClanLib_Core_LIBRARY} ${ClanLib_App_LIBRARY} ${ClanLib_Sound_LIBRARY} ${ClanLib_Display_LIBRARY} ${ClanLib_GL_LIBRARY}  ${ClanLib_GUI_LIBRARY}  -lclan23CSSLayout -lGL -lX11 -lfreetype -lfontconfig -pthread OpenCL)
#TARGET_LINK_LIBRARIES(clangl ${ClanLib_Core_LIBRARY} ${ClanLib_App_LIBRARY} ${ClanLib_Sound_LIBRARY} ${ClanLib_Display_LIBRARY} ${ClanLib_GL_LIBRARY}  ${ClanLib_GUI_LIBRARY}  -lclan23CSSLayout -lGL -lX11 -lfreetype -lfontconfig -pthread)
#TARGET_LINK_LIBRARIES(shooter2 ${ClanLib_Core_LIBRARY} ${ClanLib_App_LIBRARY} ${ClanLib_Sound_LIBRARY} ${ClanLib_Display_LIBRARY} ${ClanLib_GL_LIBRARY} -lX11 -lfreetype -lfontconfig -pthread)
#TARGET_LINK_LIBRARIES(rad ${ClanLib_Core_LIBRARY} ${ClanLib_App_LIBRARY} ${ClanLib_Sound_LIBRARY} ${ClanLib_Display_LIBRARY} ${ClanLib_GL_LIBRARY}  -lX11 -lfreetype -lfontconfig -pthread)
