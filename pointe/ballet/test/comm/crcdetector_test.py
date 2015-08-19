import unittest
import numpy as np
import ballet

class CRCDetectorTest(unittest.TestCase):

    def notatest(self):
        print "not run"

    def test_constructor(self):

        # default CRC generator polynomial
        polynomial = np.array(np.matrix('1,0,0,0,1,'
            '0,0,0,0,0,0,1,0,0,0,0,1'))
        
        # default constructor
        crcDetector = ballet.comm.CRCDetector()

        # ensure object has default values
        self.assertTrue(np.all(crcDetector.Polynomial == polynomial))
        self.assertFalse(crcDetector.ReflectChecksums)

        # CRC generator polynomial 0x3D65
        polynomial = np.array(np.matrix('1,0,0,1,1,1,'
            '1,0,1,0,1,1,0,0,1,0,1'))

        # constructor with Polynomial property
        crcDetector = ballet.comm.CRCDetector('0x3D65')

        # check object has correct property values
        self.assertTrue(np.all(crcDetector.Polynomial == polynomial))
        self.assertFalse(crcDetector.ReflectChecksums)

        # set ReflectChecksums property
        crcDetector.ReflectChecksums = True
        
        # check object has correct property values
        self.assertTrue(crcDetector.ReflectChecksums)

    def test_getProperty(self):

        # construct default CRCDetector object
        crcDetector = ballet.comm.CRCDetector()

        # get Polynomial property
        self.assertTrue(isinstance(crcDetector.Polynomial,np.ndarray))
        self.assertTrue(isinstance(crcDetector.ReflectChecksums,bool))

    def test_setProperty(self):

        # crc generator polynomial
        polynomial = np.array(np.matrix('1,0,0,1'))

        # construct CRCDetector object
        crcDetector = ballet.comm.CRCDetector()

        # set properties
        crcDetector.Polynomial = polynomial
        crcDetector.ReflectChecksums = True
        
        # assert
        self.assertTrue(np.all(crcDetector.Polynomial==polynomial))
        self.assertTrue(crcDetector.ReflectChecksums)

    def test_isLocked(self):

        # construct default CRCDetector object
        crcDetector = ballet.comm.CRCDetector()

        # check system object unlocked
        self.assertFalse(crcDetector.isLocked())

        # random binary input array
        X = np.random.randint(2,size=(100,1))
        err = crcDetector.step(X)

        # ensure CRC failure
        self.assertTrue(crcDetector.isLocked())

    def test_step(self):
        
        # construct default CRCDetector object
        crcDetector = ballet.comm.CRCDetector()

        # random binary input array
        X = np.random.randint(2,size=(100,1))
        err = crcDetector.step(X)

        # ensure CRC failure
        self.assertTrue(err)
