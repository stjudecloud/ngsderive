# Overview

The `instrument` subcommand will attempt to backward compute the machine that
generated a NGS file using (1) the instrument id(s) and (2) the flowcell id(s). 

# Limitations

* This command may not comprehensively detect the correct machines as there is
  no published catalog of Illumina serial numbers. As we encounter more serial
  numbers in practice, we update this code.
