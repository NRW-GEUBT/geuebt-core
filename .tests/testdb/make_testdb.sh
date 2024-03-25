#!/usr/bin/env bash

SCRIPT_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)

DB="geuebt-test"

mongoimport --db $DB --jsonArray --file $SCRIPT_PATH/data/runs.json
mongoimport --db $DB --jsonArray --file $SCRIPT_PATH/data/isolates.json
mongoimport --db $DB --jsonArray --file $SCRIPT_PATH/data/clusters.json