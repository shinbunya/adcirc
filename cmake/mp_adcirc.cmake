if(BUILD_MP_ADCIRC)

set(MP_ADCIRC_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sizes.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/constants.F
      ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/boundaries.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/global_3dvs.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/messenger.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/vew1d.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/harm.F
      ${CMAKE_CURRENT_SOURCE_DIR}/wind/vortex.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wind.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/hashtable.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owiwind.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owiwind_netcdf.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/rs2.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owi_ice.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/itpackv.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/nodalattr.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/globalio.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/subdomain.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/gwce.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wetdry.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/momentum.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdfio.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/control.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/xdmfio.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/writer.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/write_output.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/couple2swan.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/adcirc.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/weir_boundary.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/read_input.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/cstart.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/hstart.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/timestep.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/vsmy.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/transport.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/driver.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/sponge_layer.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/quadrature.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/couple2baroclinic3D.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/mp_messenger.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/mp_interface.F)

  add_executable(mp_adcirc ${MP_ADCIRC_SOURCES})

  addcompilerflags(mp_adcirc)
  addmpmd(mp_adcirc)
  addlibmkdir(mp_adcirc)
  addlibversion(mp_adcirc)

  add_dependencies(mp_adcirc version mkdir)

  install(TARGETS mp_adcirc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_MP_ADCIRC)
