add_subdirectory(math)
add_subdirectory(outputProcess)
add_subdirectory(berny)
add_subdirectory(chemicalUnits)




pybind11_add_module(MoleKing_Utils main.cpp)

target_link_libraries(MoleKing_Utils 
PUBLIC math
PUBLIC chemicalUnits
PUBLIC outputProcess
PUBLIC berny)