#!/bin/bash
set -euo pipefail
IFS=$'\n\t'
trap 'echo FATAL ERROR EXIT CODE $? AT $0:$LINENO' ERR

sudo docker build --progress plain --file Dockerfile -t shcsgwas:latest .

# Run the GWAS pipeline
sudo docker run \
    --volume "./output:/root/output" \
    --volume "./input:/root/input" \
    --volume "./src:/root/src" \
    --volume "./utility:/root/utility" \
    shcsgwas:latest /bin/bash /root/src/entry_script.sh

# Run the heritability pipeline
sudo docker run \
    --volume "./output:/root/output" \
    --volume "./input:/root/input" \
    --volume "./src:/root/src" \
    --volume "./utility:/root/utility" \
    shcsgwas:latest /bin/bash /root/src/entry_script_heri.sh
    