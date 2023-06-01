# CMakeLists.txt for KDetSim
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/config
)

include(ParseVersion)

project(KDetSim VERSION ${KDETSIM_VERSION})

INCLUDE(GNUInstallDirs)
include(CompilerConfig)
include(ROOT)

#---Define includes
include_directories(inc)
set(Project_Includes
    K3D.h
    KDetector.h
    KField.h
    KGeometry.h
    KImplant2D.h
    KImplant3D.h
    KMaterial.h
    KMesh.h
    KPad.h
    KPixel.h
    KStrip.h
    KStruct.h

    fizika.h
    nrutil.h
)

#---Define sources
set(Project_Sources
    src/K3D.cxx
    src/KDetector.cxx
    src/KField.cxx
    src/KGeometry.cxx
    src/KImplant2D.cxx
    src/KImplant3D.cxx
    src/KMaterial.cxx
    src/KMesh.cxx
    src/KMinuit.cxx
    src/KPad.cxx
    src/KPixel.cxx
    src/KStrip.cxx
    src/KStruct.cxx
    src/nrutil.cxx
)

#---ROOT requirements
ROOT_GENERATE_DICTIONARY(
    G__${PROJECT_NAME}

    ${Project_Includes}

    MODULE ${PROJECT_NAME}
    LINKDEF inc/KDetSim_LinkDef.h
)
set(Project_Sources ${Project_Sources} G__${PROJECT_NAME}.cxx)

#---Create a shared library with geneated dictionary
add_library(${PROJECT_NAME} SHARED
    ${Project_Sources}
)

list(TRANSFORM Project_Includes PREPEND "inc/")

set_target_properties(${PROJECT_NAME} PROPERTIES
    PUBLIC_HEADER "${Project_Includes}")

target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES})

#---Install
install(TARGETS ${PROJECT_NAME} EXPORT KDetSimTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)
install(FILES
    ${CMAKE_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm
    ${CMAKE_BINARY_DIR}/lib${PROJECT_NAME}.rootmap
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

include(Package)
