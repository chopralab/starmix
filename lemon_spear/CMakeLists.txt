add_executable(complete_score_pdb complete_score_pdb.cpp)

target_link_libraries(complete_score_pdb PRIVATE
    lemon
    chemfiles
    spear
    pthread
)

set_target_properties(complete_score_pdb PROPERTIES INSTALL_RPATH "\$ORIGIN/../lib")

install(TARGETS complete_score_pdb RUNTIME DESTINATION bin)


add_executable(idatm_protein_name idatm_protein_name.cpp)

target_link_libraries(idatm_protein_name PRIVATE
    lemon
    chemfiles
    spear
    pthread
)

set_target_properties(idatm_protein_name PROPERTIES INSTALL_RPATH "\$ORIGIN/../lib")

install(TARGETS idatm_protein_name RUNTIME DESTINATION bin)

add_executable(idatm_idatm idatm_idatm.cpp)

target_link_libraries(idatm_idatm PRIVATE
    lemon
    chemfiles
    spear
    pthread
)

set_target_properties(idatm_idatm PROPERTIES INSTALL_RPATH "\$ORIGIN/../lib")

install(TARGETS idatm_idatm RUNTIME DESTINATION bin)
