import csv
import itertools
import pysam
import re
import sys

import logging
from collections import defaultdict

from ..utils import NGSFile

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

instrument_ids = {
    "^HWI-M[0-9]{4}$": ["MiSeq"],
    "^HWUSI": ["Genome Analyzer IIx"],
    "^M[0-9]{5}$": ["MiSeq"],
    "^HWI-C[0-9]{5}$": ["HiSeq 1500"],
    "^C[0-9]{5}$": ["HiSeq 1500"],
    "^HWI-ST[0-9]{3,5}(_[0-9]{9})?$": ["HiSeq 2000"],
    "^HWI-D[0-9]{5}$": ["HiSeq 2000", "HiSeq 2500"],
    "^A[0-9]{5}$": ["NovaSeq"],
    "^D[0-9]{5}$": ["HiSeq 2500"],
    "^J[0-9]{5}$": ["HiSeq 3000"],
    "^K[0-9]{5}$": ["HiSeq 3000", "HiSeq 4000"],
    "^E[0-9]{5}$": ["HiSeq X"],
    "^N[0-9]{5}$": ["NextSeq"],
    "^NB[0-9]{6}$": ["NextSeq"],
    "^NS[0-9]{6}$": ["NextSeq"],
    "^MN[0-9]{5}$": ["MiniSeq"],
}

flowcell_ids = {
    "^C[A-Z0-9]{4}ANXX$": [
        "HiSeq 1500",
        "HiSeq 2000",
        "HiSeq 2500",
    ],  # High Output (8-lane) v4 flow cell
    "^C[A-Z0-9]{4}ACXX$": [
        "HiSeq 1000",
        "HiSeq 1500",
        "HiSeq 2000",
        "HiSeq 2500",
    ],  # High Output (8-lane) v3 flow cell
    "^D[A-Z0-9]{4}ACXX$": [
        "HiSeq 1000",
        "HiSeq 1500",
        "HiSeq 2000",
        "HiSeq 2500",
    ],  # High Output (8-lane) v3 flow cell
    "^H[A-Z0-9]{4}ADXX$": [
        "HiSeq 1500",
        "HiSeq 2000",
        "HiSeq 2500",
    ],  # Rapid Run (2-lane) v1 flow cell
    "^H[A-Z0-9]{4}BCXX$": [
        "HiSeq 1500",
        "HiSeq 2500",
    ],  # Rapid Run (2-lane) v2 flow cell
    "^H[A-Z0-9]{4}BCXY$": [
        "HiSeq 1500",
        "HiSeq 2500",
    ],  # Rapid Run (2-lane) v2 flow cell
    "^H[A-Z0-9]{4}BBXX$": ["HiSeq 4000"],  # (8-lane) v1 flow cell
    "^H[A-Z0-9]{4}BBXY$": ["HiSeq 4000"],  # (8-lane) v1 flow cell
    "^H[A-Z0-9]{4}CCXX$": ["HiSeq X"],  # (8-lane) flow cell
    "^H[A-Z0-9]{4}CCXY$": ["HiSeq X"],  # (8-lane) flow cell
    "^H[A-Z0-9]{4}ALXX$": ["HiSeq X"],  # (8-lane) flow cell
    "^H[A-Z0-9]{4}BGX[A-Z,0-9]$": ["NextSeq"],  # High output flow cell
    "^H[A-Z0-9]{4}AFXX$": ["NextSeq"],  # Mid output flow cell
    "^H[A-Z0-9]{5}RXX$": ["NovaSeq"],  # S1 flow cell
    "^H[A-Z0-9]{5}RXX$": ["NovaSeq"],  # SP flow cell
    "^H[A-Z0-9]{5}MXX$": ["NovaSeq"],  # S2 flow cell
    "^H[A-Z0-9]{5}SXX$": ["NovaSeq"],  # S4 flow cell
    "^A[A-Z0-9]{4}$": ["MiSeq"],  # MiSeq flow cell
    "^B[A-Z0-9]{4}$": ["MiSeq"],  # MiSeq flow cell
    "^D[A-Z0-9]{4}$": ["MiSeq"],  # MiSeq nano flow cell
    # "^D[A-Z0-9]{4}$" : ["HiSeq 2000", "HiSeq 2500"],                                 # Unknown HiSeq flow cell examined in SJ data
    "^G[A-Z0-9]{4}$": ["MiSeq"],  # MiSeq micro flow cell
}

upgrade_sets = [(set(["HiSeq 2000", "HiSeq 2500"]), ["HiSeq 2000", "HiSeq 2500"])]


def derive_instrument_from_iid(iid):
    matching_instruments = set()

    for pattern in instrument_ids.keys():
        if re.search(pattern, iid):
            matching_instruments |= set(instrument_ids[pattern])

    return matching_instruments


def predict_instrument_from_iids(iids):
    possible_instruments_by_iid = set()
    detected_at_least_one_instrument = False

    for iid in iids:
        if len(possible_instruments_by_iid) == 0:
            possible_instruments_by_iid = derive_instrument_from_iid(iid)
            detected_at_least_one_instrument = True
        else:
            possible_instruments_by_iid &= derive_instrument_from_iid(iid)

    return possible_instruments_by_iid, detected_at_least_one_instrument


def derive_instrument_from_fcid(fcid):
    matching_instruments = set()
    detected_at_least_one_instrument = False

    for pattern in flowcell_ids.keys():
        if re.search(pattern, fcid):
            matching_instruments |= set(flowcell_ids[pattern])

    return matching_instruments


def predict_instrument_from_fcids(fcids):
    possible_instruments_by_fcid = set()
    detected_at_least_one_instrument = False

    for fcid in fcids:
        if len(possible_instruments_by_fcid) == 0:
            possible_instruments_by_fcid = derive_instrument_from_fcid(fcid)
            detected_at_least_one_instrument = True
        else:
            possible_instruments_by_fcid &= derive_instrument_from_fcid(fcid)

    return possible_instruments_by_fcid, detected_at_least_one_instrument


def resolve_instrument(
    possible_instruments_by_iid,
    possible_instruments_by_fcid,
    at_least_one_instrument_detected,
    malformed_read_names_detected,
):
    if len(possible_instruments_by_iid) == 0 and len(possible_instruments_by_fcid) == 0:
        if at_least_one_instrument_detected:
            return (
                set(["multiple instruments"]),
                "unknown confidence",
                "multiple instruments were detected in this sample",
            )
        elif malformed_read_names_detected:
            return (
                set(["unknown"]),
                "no confidence",
                "read names not in Illumina format",
            )
        else:
            return set(["unknown"]), "no confidence", "no match"

    if len(possible_instruments_by_iid) == 0:
        confidence = "medium confidence"
        if len(possible_instruments_by_fcid) > 1:
            confidence = "low confidence"
        return possible_instruments_by_fcid, confidence, "flowcell id"

    if len(possible_instruments_by_fcid) == 0:
        confidence = "medium confidence"
        if len(possible_instruments_by_fcid) > 1:
            confidence = "low confidence"
        return possible_instruments_by_iid, confidence, "instrument id"

    overlapping_instruments = possible_instruments_by_iid & possible_instruments_by_fcid
    if len(overlapping_instruments) == 0:
        return (
            set(["conflicting evidence"]),
            "high confidence",
            "Case needs triaging: {} by iid, {} by fcid".format(
                " or ".join(possible_instruments_by_iid),
                " or ".join(possible_instruments_by_fcid),
            ),
        )
    else:
        return (
            overlapping_instruments,
            "high confidence",
            "instrument id and flowcell id",
        )


def main(ngsfiles, outfile=sys.stdout, delimiter="\t", n_samples=10000):
    writer = csv.DictWriter(
        outfile,
        fieldnames=["File", "Instrument", "Confidence", "Basis"],
        delimiter=delimiter,
    )
    writer.writeheader()
    outfile.flush()

    if n_samples < 1:
        n_samples = None

    for ngsfilepath in ngsfiles:
        try:
            ngsfile = NGSFile(ngsfilepath)
        except FileNotFoundError:
            result = {
                "File": ngsfilepath,
                "Instrument": "File not found.",
                "Confidence": "N/A",
                "Basis": "N/A",
            }

            writer.writerow(result)
            outfile.flush()
            continue

        instruments = set()
        flowcells = set()
        malformed_read_names = False

        # accumulate instrument and flowcell IDs
        for read in itertools.islice(ngsfile, n_samples):
            parts = read["query_name"].split(":")
            if len(parts) != 7:  # not Illumina format
                malformed_read_names = True
                continue
            iid, fcid = parts[0], parts[2]
            instruments.add(iid)
            flowcells.add(fcid)

        (
            possible_instruments_by_iid,
            detected_instrument_by_iid,
        ) = predict_instrument_from_iids(instruments)
        (
            possible_instruments_by_fcid,
            detected_instrument_by_fcid,
        ) = predict_instrument_from_fcids(flowcells)

        instruments, confidence, based_on = resolve_instrument(
            possible_instruments_by_iid,
            possible_instruments_by_fcid,
            detected_instrument_by_iid | detected_instrument_by_fcid,
            malformed_read_names,
        )
        for upgrade_set in upgrade_sets:
            if instruments.issubset(upgrade_set[0]):
                instruments = upgrade_set[1]
                break

        result = {
            "File": ngsfilepath,
            "Instrument": " or ".join(instruments),
            "Confidence": confidence,
            "Basis": based_on,
        }

        writer.writerow(result)
        outfile.flush()
