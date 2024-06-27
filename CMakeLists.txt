cmake_minimum_required(VERSION 3.1)

set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)

set(PROJECT_NAME FSMSequenceOne)
set(PROJECT_DESCRIPTION "FSMSequenceOne")
set(PROJECT_URL "")

project(${PROJECT_NAME} LANGUAGES CXX VERSION 1.0.0)

find_package(mc_rtc REQUIRED)
find_package(mc_iam REQUIRED)

set(SHARE_DESTINATION "${CMAKE_INSTALL_PREFIX}/share/toss_scene")
set(MESHES_DESTINATION "${SHARE_DESTINATION}/meshes")
set(MODULES_DESTINATION "${SHARE_DESTINATION}/modules")
set(URDF_DESTINATION "${SHARE_DESTINATION}/urdf")

macro(install_aliase ROBOT)
  set(ALIASE_OUT "${CMAKE_CURRENT_BINARY_DIR}/aliases/${ROBOT}.yaml")
  configure_file(aliases/${ROBOT}.in.yaml "${ALIASE_OUT}")
  install(FILES "${ALIASE_OUT}" DESTINATION "${MC_ROBOTS_ALIASES_DIRECTORY}")
endmacro()

# VANDERLANDE WORKBENCH INSTALLATION
set(URDF_WORKBENCH "${CMAKE_CURRENT_BINARY_DIR}/urdf/workbench.urdf")
configure_file(urdf/workbench.urdf "${URDF_WORKBENCH}" @ONLY)
install(FILES "${URDF_WORKBENCH}" DESTINATION "${URDF_DESTINATION}")
install_aliase(workbench)

# VANDERLANDE CONVEYOR
set(URDF_CONVEYOR002 "${CMAKE_CURRENT_BINARY_DIR}/urdf/conveyor002.urdf")
configure_file(urdf/conveyor002.urdf "${URDF_CONVEYOR002}" @ONLY)
install(FILES "${URDF_CONVEYOR002}" DESTINATION "${URDF_DESTINATION}")
install_aliase(conveyor002)

# VANDERLANDE TOTE
set(URDF_TOTE "${CMAKE_CURRENT_BINARY_DIR}/urdf/tote.urdf")
configure_file(urdf/tote.urdf "${URDF_TOTE}" @ONLY)
install(FILES "${URDF_TOTE}" DESTINATION "${URDF_DESTINATION}")
install_aliase(tote)

# Install surfaces description
install(DIRECTORY rsdf DESTINATION "${SHARE_DESTINATION}" FILES_MATCHING PATTERN "*.rsdf")

# Install convex
install(DIRECTORY convex DESTINATION "${SHARE_DESTINATION}" FILES_MATCHING PATTERN "*-ch.txt")

# Install meshes
install(DIRECTORY meshes DESTINATION "${SHARE_DESTINATION}")

# BOXES INSTALLATION
install(FILES urdf/box1.urdf DESTINATION "${URDF_DESTINATION}")
install(FILES urdf/box2.urdf DESTINATION "${URDF_DESTINATION}")
install(FILES urdf/box3.urdf DESTINATION "${URDF_DESTINATION}")
install(FILES urdf/box4.urdf DESTINATION "${URDF_DESTINATION}")
install(FILES urdf/box5.urdf DESTINATION "${URDF_DESTINATION}")
install(FILES urdf/box013.urdf DESTINATION "${URDF_DESTINATION}")

install_aliase(toss_box1)
install_aliase(toss_box2)
install_aliase(toss_box3)
install_aliase(toss_box4)
install_aliase(toss_box5)
install_aliase(box013)

add_subdirectory(src)

set(AROBASE "@")
configure_file(etc/FSMSequenceOne.in.yaml "${CMAKE_CURRENT_BINARY_DIR}/etc/FSMSequenceOne.yaml")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/etc/FSMSequenceOne.yaml" DESTINATION "${MC_CONTROLLER_RUNTIME_INSTALL_PREFIX}/etc")

