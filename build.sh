#!/usr/bin/env bash
set -euo pipefail

# Deprecated compatibility wrapper. Prefer scripts/dev-build.sh.
exec "$(dirname "$0")/scripts/dev-build.sh" --serial --links "$@"
