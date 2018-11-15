#
# Barcode
# Copyright E.G.P. Bos and F.S. Kitaura
#
# Distributed under the terms of the MIT License.
# The full license is in the file LICENSE, distributed with this software.
#

add_executable(barcode main.cc)

target_link_libraries(barcode PRIVATE barlib)
TARGET_ADD_FFTW(barcode)
if (MULTITHREAD)
    target_link_libraries(barcode PRIVATE OpenMP::OpenMP_CXX)
endif()

# make a link for the most recent build
# if Mercurial is present, use it to add the revision number to the output executable
# N.B.: this causes the barcode binary output file to not be present, so every build will always do a new linking step!
#find_package(Hg)
#if(HG_FOUND)
#    get_filename_component(OUTPUT_FN $<TARGET_FILE_NAME:barcode> NAME)
#    add_custom_command(TARGET barcode
#            POST_BUILD
#            COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/hg_rev_filename.sh ${PROJECT_BINARY_DIR} ${OUTPUT_FN}
#            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#            COMMENT "Renamed output binary to include Mercurial revision information and created symbolic link barcode_latest_build to new path.")
#endif(HG_FOUND)

#--------------------------------------- installation
install(TARGETS barcode
        RUNTIME DESTINATION bin)
