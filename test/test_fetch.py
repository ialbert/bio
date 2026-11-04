import io
import zipfile
from unittest.mock import patch

from biorun import fetch


def make_dataset():
    data = io.BytesIO()
    root = "ncbi_dataset/data/GCF_000005845.2"
    with zipfile.ZipFile(data, "w") as archive:
        archive.writestr(f"{root}/GCF_000005845.2_ASM584v2_genomic.fna", ">seq\nACGT\n")
        archive.writestr(f"{root}/genomic.gbff", "LOCUS       seq\n")
        archive.writestr(f"{root}/genomic.gff", "##gff-version 3\n")
        archive.writestr(f"{root}/genomic.gtf", "# gtf\n")
        archive.writestr(f"{root}/cds_from_genomic.fna", ">cds\nACG\n")
        archive.writestr(f"{root}/rna.fna", ">rna\nACG\n")
        archive.writestr(f"{root}/protein.faa", ">protein\nM\n")
        archive.writestr(f"{root}/sequence_report.jsonl", "{}\n")
    return data.getvalue()


class Response:
    def __init__(self, data):
        self.data = data
        self.headers = {"content-length": str(len(data))}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None

    def raise_for_status(self):
        return None

    def iter_content(self, chunk_size):
        yield self.data


def test_is_gcf():
    assert fetch.is_GCF("GCF_000005845.2")
    assert fetch.is_GCF("GCA_000001405")
    assert not fetch.is_GCF("GCF_bad")
    assert not fetch.is_GCF("prefix_GCF_000005845.2")
    assert fetch.format_mb(1234 * 1024 * 1024) == "1,234Mb"
    assert fetch.assembly_dir("genome.fna") == "fasta"
    assert fetch.assembly_dir("protein.faa") == "fasta"
    assert fetch.assembly_dir("ann.gff") == "gff"
    assert fetch.assembly_dir("ann.gtf") == "gff"


@patch("biorun.fetch.requests.get")
def test_fetch_gcf(mock_get, tmp_path, capsys, monkeypatch):
    mock_get.return_value = Response(make_dataset())
    monkeypatch.chdir(tmp_path)

    fetch.fetch_gcf(["gcf_000005845.2"])

    expected = {
        tmp_path / "fasta/GCF_000005845.2.genome.fna",
        tmp_path / "gff/GCF_000005845.2.ann.gff",
        tmp_path / "gff/GCF_000005845.2.ann.gtf",
    }
    assert all(path.exists() for path in expected)
    assert (tmp_path / "fasta/GCF_000005845.2.genome.fna").read_text() == ">seq\nACGT\n"
    output = {f"# {path.relative_to(tmp_path)} (0Mb)" for path in expected}
    assert set(capsys.readouterr().out.splitlines()) == output

    _, kwargs = mock_get.call_args
    assert kwargs["headers"] == {"Accept": "application/zip"}
    assert len(kwargs["params"]) == 3


@patch("biorun.fetch.requests.get")
def test_fetch_gcf_rename(mock_get, tmp_path, capsys):
    mock_get.return_value = Response(make_dataset())

    fetch.fetch_gcf(["GCF_000005845.2"], include="gff,cds", out=tmp_path, rename="e-coli")

    gff = tmp_path / "e-coli.ann.gff"
    cds = tmp_path / "e-coli.cds.fna"
    assert {path.name for path in tmp_path.iterdir()} == {gff.name, cds.name}
    output = {f"# {gff} (0Mb)", f"# {cds} (0Mb)"}
    assert set(capsys.readouterr().out.splitlines()) == output
