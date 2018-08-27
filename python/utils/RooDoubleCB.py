import os
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

r.gSystem.Load(os.path.join(os.environ['CHIB_CHIC_POLFW_DIR'],
                            'general', 'RooDoubleCB', 'RooDoubleCB.so'))
