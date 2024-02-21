#!/bin/bash
ARGS=$@

eval "$(/usr/local/bin/brew shellenv)"
arch -x86_64 /bin/bash -c "${ARGS}"
