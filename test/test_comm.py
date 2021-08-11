from biorun.scripts import comm
import biorun.libs.placlib as plac
import os

join = os.path.join

# The path to the current file.
CURR_DIR = os.path.dirname(__file__)

# The default data directory.
DATA_DIR = join(CURR_DIR, "data")

FNAME1= join(DATA_DIR, "file1.txt")
FNAME2= join(DATA_DIR, "file2.txt")

def test_comm_1(capsys):
    params = dict(
        file1=FNAME1,
        file2=FNAME2,
    )
    comm.main(**params)
    stream = capsys.readouterr().out
    assert stream == 'A\nC\n'

def test_comm_2(capsys):
    params = dict(
        file1=FNAME1,
        file2=FNAME2,
        uniq1=True,
    )
    comm.main(**params)
    stream = capsys.readouterr().out
    assert stream == 'B\n'

def test_comm_3(capsys):
    params = dict(
        file1=FNAME1,
        file2=FNAME2,
        uniq2=True
    )
    comm.main(**params)
    stream = capsys.readouterr().out
    assert stream == 'D\n'

def test_comm_4(capsys):
    params = dict(
        file1=FNAME1,
        file2=FNAME2,
        union=True
    )
    comm.main(**params)
    stream = capsys.readouterr().out
    assert stream == 'A\nC\nB\nD\n'


