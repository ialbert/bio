from biorun.methods import comm
import biorun.libs.placlib as plac
import os


FNAME1="test/data/file1.txt"
FNAME2="test/data/file2.txt"


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


