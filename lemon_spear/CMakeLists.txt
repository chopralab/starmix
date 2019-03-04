add_executable(complete_score_pdb complete_score_pdb.cpp)

target_link_libraries(complete_score_pdb PRIVATE
    lemon
    chemfiles
    spear
    pthread
)

set_target_properties(complete_score_pdb PROPERTIES INSTALL_RPATH "\$ORIGIN/../lib")

install(TARGETS complete_score_pdb RUNTIME DESTINATION bin)
