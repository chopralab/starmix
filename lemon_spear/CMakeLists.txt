function(add_cpp_prog _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    add_executable(${_name_} ${_file_})

    target_link_libraries(${_name_} PRIVATE
        lemon
        spear
    )

    set_target_properties(${_name_} PROPERTIES INSTALL_RPATH "\$ORIGIN/../lib")
    install(TARGETS ${_name_} RUNTIME DESTINATION bin)
endfunction()

add_cpp_prog(complete_score_pdb.cpp)
add_cpp_prog(idatm_protein_name.cpp)
add_cpp_prog(idatm_idatm.cpp)
