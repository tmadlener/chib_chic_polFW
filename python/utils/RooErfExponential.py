import os
import logging
logging.basicConfig(level=logging.ERROR,
                    format='%(levelname)s - %(funcName)s: %(message)s')

from utils.root_lib_loading import try_load_lib

lib_path = os.path.join(os.environ['CHIB_CHIC_POLFW_DIR'],
                        'general', 'RooErfExponential', 'RooErfExponential.so')

try_load_lib(lib_path, 'RooErfExponential')
