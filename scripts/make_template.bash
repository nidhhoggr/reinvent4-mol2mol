#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_SCRIPT_NAME=${1:-libeinvent_rl}
python ${SCRIPT_DIR}/make_template.py \
    ${SCRIPT_DIR}/../configs/${CONFIG_SCRIPT_NAME}.toml \
    ${SCRIPT_DIR}/../configs/${CONFIG_SCRIPT_NAME}.template.toml
