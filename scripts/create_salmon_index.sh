#!/bin/bash

cd /projects/mcgaugheyd/genome/GRCh38

~/bin/Salmon-0.7.2_linux_x86_64/bin/./salmon index --transcripts gencode.v25.pc_transcripts.fa.gz --index salmon_0.7.2_gencode_v25 --gencode --type quasi --perfectHash -p 16
