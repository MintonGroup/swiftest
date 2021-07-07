#!/bin/bash
gprof ./swifter_rmvs | /home/daminton/git/gprof2dot/gprof2dot.py | dot -Tpng -o swifter_profile.png
