#!/usr/bin/env bash

DB="geuebt-test"

mongosh --eval "use $DB" --eval  "db.dropDatabase()"