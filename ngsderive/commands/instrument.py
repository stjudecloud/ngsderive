import csv
import itertools
import logging
import re

from ..utils import NGSFile

logger = logging.getLogger("instrument")

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
    "^H[A-Z0-9]{5}RXX$": ["NovaSeq"],  # S1/SP flow cell
    "^H[A-Z0-9]{5}MXX$": ["NovaSeq"],  # S2 flow cell
    "^H[A-Z0-9]{5}SXX$": ["NovaSeq"],  # S4 flow cell
    "^A[A-Z0-9]{4}$": ["MiSeq"],  # MiSeq flow cell
    "^B[A-Z0-9]{4}$": ["MiSeq"],  # MiSeq flow cell
    "^D[A-Z0-9]{4}$": ["MiSeq"],  # MiSeq nano flow cell
    # "^D[A-Z0-9]{4}$" : ["HiSeq 2000", "HiSeq 2500"],  # Unknown HiSeq flow cell examined in SJ data
    "^G[A-Z0-9]{4}$": ["MiSeq"],  # MiSeq micro flow cell
}

upgrade_sets = [(set(["HiSeq 2000", "HiSeq 2500"]), ["HiSeq 2000", "HiSeq 2500"])]


def derive_instrument_from_iid(iid):
    matching_instruments = set()

    for pattern, instruments in instrument_ids.items():
        if re.search(pattern, iid):
            matching_instruments |= set(instruments)

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

    for pattern, instruments in flowcell_ids.items():
        if re.search(pattern, fcid):
            matching_instruments |= set(instruments)

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
        if malformed_read_names_detected:
            return (
                set(["unknown"]),
                "no confidence",
                "read names not in Illumina format",
            )
        return set(["unknown"]), "no confidence", "no match"

    if len(possible_instruments_by_iid) == 0:
        if not malformed_read_names_detected:
            confidence = "medium confidence"
            if len(possible_instruments_by_fcid) > 1:
                confidence = "low confidence"
            return possible_instruments_by_fcid, confidence, "flowcell id"
        confidence = "low confidence"
        return possible_instruments_by_fcid, confidence, "flowcell id"

    if len(possible_instruments_by_fcid) == 0:
        if not malformed_read_names_detected:
            confidence = "medium confidence"
            if len(possible_instruments_by_iid) > 1:
                confidence = "low confidence"
            return possible_instruments_by_iid, confidence, "instrument id"
        confidence = "low confidence"
        return (
            possible_instruments_by_iid,
            confidence,
            "instrument id",
        )

    overlapping_instruments = possible_instruments_by_iid & possible_instruments_by_fcid
    if len(overlapping_instruments) == 0:
        return (
            set(["conflicting evidence"]),
            "high confidence",
            f"Case needs triaging: {' or '.join(possible_instruments_by_iid)} by iid, {' or '.join(possible_instruments_by_fcid)} by fcid",
        )
    return (
        overlapping_instruments,
        "high confidence",
        "instrument id and flowcell id",
    )


def main(ngsfiles, outfile, n_reads):
    writer = csv.DictWriter(
        outfile,
        fieldnames=["File", "Instrument", "Confidence", "Basis"],
        delimiter="\t",
    )
    writer.writeheader()
    outfile.flush()

    if n_reads < 1:
        n_reads = None

    for ngsfilepath in ngsfiles:
        try:
            ngsfile = NGSFile(ngsfilepath)
        except FileNotFoundError:
            result = {
                "File": ngsfilepath,
                "Instrument": "Error opening file.",
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
        try:
            for read in itertools.islice(ngsfile, n_reads):
                parts = read["query_name"].split(":")
                if len(parts) != 7:  # not Illumina format
                    malformed_read_names = True
                    iid = parts[0]  # attempt to recover machine name
                    instruments.add(iid)
                    for rg in ngsfile.handle.header.to_dict()["RG"]:
                        if rg["ID"] == read["read_group"]:
                            if "PU" in rg:
                                flowcells.add(rg["PU"])
                            if "PM" in rg:
                                instruments.add(rg["PM"])
                    continue
                iid, fcid = parts[0], parts[2]
                instruments.add(iid)
                flowcells.add(fcid)
        except KeyError:  # no RG tag is present
            result = {
                "File": ngsfilepath,
                "Instrument": "unknown",
                "Confidence": "no confidence",
                "Basis": "no RG tag present",
            }
            writer.writerow(result)
            outfile.flush()
            continue

        if malformed_read_names:
            logger.warning(
                "Encountered read names not in Illumina format. Recovery attempted."
            )
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
