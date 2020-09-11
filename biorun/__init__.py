
import os


__CURRENT_DIR = os.path.dirname(__file__)
__DUMP = os.path.join(__CURRENT_DIR, '..', 'save_file')

VERSION = "0.0.1"

# Default directory used to store data
DUMP_DIR = os.path.abspath(os.getenv('DUMP_DIR', __DUMP))
