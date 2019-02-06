import sys
import os
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.ERROR,
                    format='%(levelname)s - %(funcName)s: %(message)s')



lib_path = os.path.join(os.environ['CHIB_CHIC_POLFW_DIR'],
                        'general', 'RooDoubleCB', 'RooDoubleCB.so')

if r.gSystem.Load(lib_path) != 0:
    # For some reason only loading the shared object below fails with an
    # undefined symbol error sometimes. This can normally be fixed by loading
    # libRooFitCore first. In principle this could be done independently of
    # first trying to load the shared object, but loading libRooFitCore
    # displays the splash text from RooFit which is something I do not like
    # entirely, so we do it only when really necessary
    logging.debug('Cannot load the RooDoubleCB shared object. Trying again after'
                  ' loading libRooFitCore first')
    r.gSystem.Load('libRooFitCore')
    if r.gSystem.Load(lib_path) != 0:
        # If we still fail to load the library then bail out. This is a bigger
        # problem than we thought
        raise ImportError('Cannot load the RooDoubleCB shared object.\nMake sure'
                          ' that it is compiled and linked against the proper '
                          'libraries.\nCheck {}/general/RooDoubleCB/ for this '
                          'and run \'make clean; make\' there and see if this '
                          'fixes things.'
                          .format(os.environ['CHIB_CHIC_POLFW_DIR']))
