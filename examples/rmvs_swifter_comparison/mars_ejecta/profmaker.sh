#!/bin/bash
gprof ./swiftest_driver | /home/daminton/git/gprof2dot/gprof2dot.py | dot -Tpng -o swiftest_profile.png
