"""
The main job runner. Register functions here.
"""
import plac
from biorun.main import run

def runner():
    plac.call(run)

if __name__ == '__main__':
    runner()
