#!/usr/bin/env bash

SCRIPT_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)

bash $SCRIPT_PATH/delete_testdb.sh
bash $SCRIPT_PATH/make_testdb.sh
