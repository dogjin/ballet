file ( GLOB BALLET_SRCS
       "src/corelib/kernel/object.cpp"
       "src/corelib/tools/random.cpp"
       "src/corelib/tools/logical.cpp"
       "src/comm/tools/bi2de.cpp"
       "src/comm/tools/de2bi.cpp"
       "src/comm/tools/oct2dec.cpp"
       "src/comm/crc/crcbase.cpp"
       "src/comm/crc/crcgenerator.cpp"
       "src/comm/crc/crcdetector.cpp"
       "src/comm/coding/convolutionalcoding.cpp"
       "src/comm/coding/convolutionalencoder.cpp"
       "src/comm/coding/viterbi.cpp"
       #"src/comm/coding/lva.cpp"
       "src/comm/coding/appdecoder.cpp"
       "src/comm/coding/turboencoder.cpp"
       "src/comm/coding/turbodecoder.cpp" )

file ( GLOB BALLET_HEADER_FILES
       "src/corelib/kernel/object.h"
       "src/corelib/global/global.h"
       "src/corelib/tools/logical.h"
       "src/corelib/tools/random.h"
       "include/ballet/ballet.h" )

file ( GLOB BALLET_COMM_HEADER_FILES
       "src/comm/tools/bi2de.h"
       "src/comm/tools/de2bi.h"
       "src/comm/tools/oct2dec.h"
       "src/comm/crc/crcdetector.h"
       "src/comm/crc/crcbase.h"
       "src/comm/crc/crcgenerator.h"
       "src/comm/coding/viterbi.h"
       "src/comm/coding/appdecoder.h"
       #"src/comm/coding/lva.h"
       "src/comm/coding/convolutionalencoder.h"
       "src/comm/coding/convolutionalcoding.h"
       "src/comm/coding/turboencoder.h"
       "src/comm/coding/turbodecoder.h"
       "include/ballet/ballet.h" )


include_directories ( ${CMAKE_SOURCE_DIR}/ballet/include
                      ${CMAKE_SOURCE_DIR}/ballet/include/ballet/1.0.0/ballet )

add_library ( ballet SHARED ${BALLET_SRCS} )

set_target_properties ( ballet PROPERTIES
    VERSION ${PACKAGE_VERSION}
    SOVERSION ${PACKAGE_VERSION_MAJOR} )

install ( TARGETS ballet
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib/static )

install ( FILES ${BALLET_HEADER_FILES}
          DESTINATION include/ballet )

install ( FILES ${BALLET_COMM_HEADER_FILES}
          DESTINATION include/ballet/comm )

add_subdirectory( test EXCLUDE_FROM_ALL )
