from CMGTools.RootTools.analyzers.DiLeptonAnalyzer import DiLeptonAnalyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.DiObject import DiElectron
from CMGTools.RootTools.physicsobjects.PhysicsObjects import Electron


class ZEleEleAnalyzer( DiLeptonAnalyzer ):

    DiObjectClass = DiElectron
    LeptonClass = Electron

    def declareHandles(self):
        super(ZEleEleAnalyzer, self).declareHandles()
        print 'ZEleEleAnalyzer.declareHandles'
        self.handles['diLeptons'] = AutoHandle(
            'cmgDiElectronSelStdLep',
            'std::vector<cmg::DiObject<cmg::Electron,cmg::Electron>>'
            )
        self.handles['leptons'] = AutoHandle(
            'cmgElectronSelStdLep',
            'std::vector<cmg::Electron>'
            )

    def testLeg1(self, leg):
        return self.testElectron(leg) and \
               super( ZEleEleAnalyzer, self).testLeg1( leg )

    def testLeg2(self, leg):
        return self.testElectron(leg) and \
               super( ZEleEleAnalyzer, self).testLeg2( leg )


