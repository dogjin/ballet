import unittest

TEST_MODULES = [
    'ballet.test.comm.crcdetector_test'
]

def all():
    return unittest.defaultTestLoader.loadTestsFromNames(TEST_MODULES)

if __name__ == '__main__':
    unittest.main(defaultTest='all')
