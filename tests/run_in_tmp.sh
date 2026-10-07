#!/bin/sh
# Run a command in a fresh scratch directory that has an output/ subdirectory,
# so the solver's files do not land in the repository. Exits with the
# command's status.
set -u
dir=$(mktemp -d)
mkdir "$dir/output"
(cd "$dir" && "$@")
status=$?
rm -rf "$dir"
exit $status
