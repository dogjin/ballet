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
                self.__d.setProperty('Polynomial',POLY)

            elif isinstance(args[0],np.ndarray):
                self.__d.setProperty('Polynomial',args[0])

            elif isinstance(args[0],CRCDetector):
                srcObj = args[0]
                self.__d.setProperty('Polynomial',srcObj.Polynomial)
                self.__d.setProperty('ReflectChecksums',srcObj.ReflectChecksums)

    def isLocked(self):
        return self.__d.isLocked()

    def step(self,X):
        return self.__d.step(X)

    def reset(self):
        return self.__d.reset()

    def release(self):
        return self.__d.release()

    def __getattr__(self,name):

        if not name in ('Polynomial','ReflectChecksums'):
            raise Exception('Unknown property %s' % name)

        return self.__d.getProperty(name)

    def __setattr__(self,name,value):

        if name in ('Polynomial','ReflectChecksums'):
            return self.__d.setProperty(name,value)

        return super(CRCDetector,self).__setattr__(name,value)
