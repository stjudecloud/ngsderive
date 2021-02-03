#!/usr/bin/env python3

import rstr
import argparse
from ngsderive.commands import instrument

instrument_ids = instrument.instrument_ids
flowcell_ids = instrument.flowcell_ids


def generate_machines(f, n=5):
    f.write("\t".join(["Instrument", "Id", "All Matching Instruments"]) + "\n")
    for key in instrument_ids.keys():
        for instrument in instrument_ids[key]:
            for i in range(n):
                iid = rstr.xeger(key)
                f.write(
                    "\t".join([instrument, iid, ",".join(instrument_ids[key])]) + "\n"
                )


def generate_flowcells(f, n=5):
    f.write("\t".join(["Flowcell", "Id", "All Matching Flowcells"]) + "\n")
    for key in flowcell_ids.keys():
        for flowcell in flowcell_ids[key]:
            for i in range(n):
                fcid = rstr.xeger(key)
                f.write("\t".join([flowcell, fcid, ",".join(flowcell_ids[key])]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Helper script to generate randomly matching IDs for instrument and flowcell patterns."
    )
    parser.add_argument(
        "--instruments-file",
        help="File to save the generated instrument ids to.",
        default="fixtures/instruments.mocked.tsv",
    )
    parser.add_argument(
        "--flowcells-file",
        help="File to save the generated flowcell ids to.",
        default="fixtures/flowcells.mocked.tsv",
    )
    args = parser.parse_args()

    with open(args.instruments_file, "w") as instruments_file, open(
        args.flowcells_file, "w"
    ) as flowcells_file:
        generate_machines(instruments_file)
        generate_flowcells(flowcells_file)
