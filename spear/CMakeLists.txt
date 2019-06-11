function(add_spear_prog _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    add_executable(${_name_} ${_file_})

    target_link_libraries(${_name_} PRIVATE
        spear
    )

    set_target_properties(${_name_} PROPERTIES INSTALL_RPATH "\$ORIGIN/../lib")
    install(TARGETS ${_name_} RUNTIME DESTINATION bin)
endfunction()

add_spear_prog(score_all_residues.cpp)
add_spear_prog(filter_carboxylic_acids.cpp)
add_spear_prog(score_poses.cpp)
add_spear_prog(score_poses_vina.cpp)
