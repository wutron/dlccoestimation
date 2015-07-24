"""

   Code for the DLC_coestimation model
   (duplications, losses, and coalescence)

"""

from __future__ import division

import sys, os

def load_deps(dirname="deps"):
    sys.path.append(os.path.realpath(
        os.path.join(os.path.dirname(__file__), dirname)))

# add pre-bundled dependencies to the python path,
# if they are not available already
try:
    import rasmus, compbio
except ImportError:
    load_deps()
    import rasmus, compbio

#=============================================================================
# constants

PROGRAM_NAME = u"DLC Coestimation"
PROGRAM_VERSION_MAJOR = 1
PROGRAM_VERSION_MINOR = 0
PROGRAM_VERSION_RELEASE = 0
PROGRAM_VERSION = (PROGRAM_VERSION_MAJOR,
                   PROGRAM_VERSION_MINOR,
                   PROGRAM_VERSION_RELEASE)

if PROGRAM_VERSION_RELEASE != 0:
    PROGRAM_VERSION_TEXT = "%d.%d.%d" % (PROGRAM_VERSION_MAJOR,
                                         PROGRAM_VERSION_MINOR,
                                         PROGRAM_VERSION_RELEASE)
else:
    PROGRAM_VERSION_TEXT = "%d.%d" % (PROGRAM_VERSION_MAJOR,
                                      PROGRAM_VERSION_MINOR)



#=============================================================================
# miscellaneous

ALIGNMENT_SEED = int("0xDEADBEEF", 16)

class NullLog (object):

    def __init__(self):
        pass

    def write(self, text):
        pass

    def flush(self):
        pass
