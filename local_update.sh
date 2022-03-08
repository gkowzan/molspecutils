#!/usr/bin/env bash
# Build updated sdist and wheel without bumping the version, upload it and
# update rotsim2d and rotsim2d_apps projects.
set -uo pipefail

die() {
    printf '%s\n' "$1" >&2
    exit 1
}

python -m build || die "build failed"
twine upload -r local dist/*

cd ../rotsim2d || die "can't find rotsim2d"
pip install --extra-index-url=http://127.0.0.1:4040/simple --no-deps --pre --upgrade --force-reinstall --no-cache-dir molspecutils
cd ../rotsim2d_apps || die "can't find rotsim2d_apps"
pip install --extra-index-url=http://127.0.0.1:4040/simple --no-deps --pre --upgrade --force-reinstall --no-cache-dir molspecutils
cd ../molspecutils
