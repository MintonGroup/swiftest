#!/bin/bash
gprof ./swifter_symba_ringmoons | /home/daminton/git/gprof2dot/gprof2dot.py | dot -Tpng -o output.png
