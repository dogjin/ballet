import unittest
import ballet

class CRCDetectorTest(unittest.TestCase):

    def notatest(self):
        print "not run"

    def test_constructor(self):
        
        # default constructor
        crcDetector = ballet.comm.CRCDetector()

        self.assertTrue(True)

        # constructor with Polynomial property
        crcDetector = ballet.comm.CRCDetector('0x1021')

    def test_getProperty(self):
        self.assertTrue(True)
