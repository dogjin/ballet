import _crcdetector

class CRCDetector(object):
    def __init__(self,*args,**kwargs):

        self.__d = _crcdetector.CRCDetector()

        nargin = len(args)

        if nargin == 1:

            if isinstance(args[0],(int,long)):
                POLY = _crcdetector.convertHexToBinary(args[0])
                self.__d.setProperty('Polynomial',POLY)

            elif isinstance(args[0],basestring):
                POLY = _crcdetector.convertHexToBinary(args[0])
                print "POLY =", POLY
                self.__d.setProperty('Polynomial',POLY)

            elif isinstance(args[0],np.ndarray):
                self.__d.setProperty('Polynomial',args[0])

            elif isinstance(args[0],CRCDetector):
                srcObj = args[0]
                self.__d.setProperty('Polynomial',srcObj.Polynomial)
                self.__d.setProperty('ReflectChecksums',srcObj.ReflectChecksums)
