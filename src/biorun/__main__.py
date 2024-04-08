"""
This module exists only to make module level execution also viable.

python -m biorun

"""

from biorun.main import router

def run():
    """
    A simple wrapper over the task_selector
    """
    router()

if __name__ == '__main__':
    run()
