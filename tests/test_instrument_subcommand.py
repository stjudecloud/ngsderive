import csv
from ngsderive.commands import instrument


def test_iid_regexes_match():
    with open("tests/fixtures/instruments.mocked.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            possible = instrument.derive_instrument_from_iid(row["Id"])
            expected = set(row["All Matching Instruments"].split(","))
            assert possible == expected


def test_fcid_regexes_match():
    with open("tests/fixtures/flowcells.mocked.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            possible = instrument.derive_instrument_from_fcid(row["Id"])
            expected = set(row["All Matching Flowcells"].split(","))
            assert possible == expected
