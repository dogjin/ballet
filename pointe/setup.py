from distutils.core import setup, Extension
import numpy

crcbase_extension = Extension('ballet.comm.crc._crcbase',
    sources = ['ballet/comm/crc/enpointe_crcbase.cpp'],
    include_dirs=[numpy.get_include()],
    libraries=['ballet','armadillo'])

crcgenerator_extension = Extension('ballet.comm.crc._crcgenerator',
    sources = ['ballet/comm/crc/enpointe_crcgenerator.cpp'],
    include_dirs=[numpy.get_include()],
    libraries=['ballet','armadillo'])

crcdetector_extension = Extension('ballet.comm.crc._crcdetector',
    sources = ['ballet/comm/crc/enpointe_crcdetector.cpp'],
    include_dirs=[numpy.get_include()],
    libraries=['ballet','armadillo'])

comm_tools_extension = Extension('ballet.comm.tools._tools',
    sources = ['ballet/comm/tools/enpointe_comm_tools.cpp'],
    include_dirs=[numpy.get_include()],
    libraries=['ballet','armadillo'])

setup (
    name = 'PointeBallet',
    version = '0.1',
    description = 'This is a signal processing library',
    packages = ['ballet','ballet.core','ballet.comm','ballet.comm.crc',
        'ballet.comm.tools','ballet.test','ballet.test.comm'],
    ext_modules = [crcbase_extension,crcdetector_extension,
        crcgenerator_extension,comm_tools_extension])
