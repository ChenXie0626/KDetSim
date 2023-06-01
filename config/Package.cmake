INCLUDE(CMakePackageConfigHelpers)

WRITE_BASIC_PACKAGE_VERSION_FILE(
    "${CMAKE_BINARY_DIR}/package/KDetSim/KDetSimConfigVersion.cmake"
    VERSION ${KDETSIM_VERSION}
    COMPATIBILITY AnyNewerVersion
)
EXPORT(EXPORT KDetSimTargets
       FILE "${CMAKE_BINARY_DIR}/package/KDetSim/KDetSimTargets.cmake"
)
CONFIGURE_FILE(config/package/KDetSimConfig.cmake
    "${CMAKE_BINARY_DIR}/package/KDetSim/KDetSimConfig.cmake"
)
SET(KDetSimConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/KDetSim)
INSTALL(EXPORT KDetSimTargets
        FILE KDetSimTargets.cmake
        DESTINATION ${KDetSimConfigPackageLocation}
)
INSTALL(
    FILES
        config/package/KDetSimConfig.cmake
        "${CMAKE_BINARY_DIR}/package/KDetSim/KDetSimConfigVersion.cmake"
    DESTINATION ${KDetSimConfigPackageLocation}
    COMPONENT Devel
)
