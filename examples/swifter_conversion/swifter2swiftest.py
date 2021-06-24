import sys
import argparse
import swiftestio as swio
"""
  Converts initial conditions files from Swifter to Swiftest
"""

if __name__ == '__main__':
   ap = argparse.ArgumentParser()
   ap.add_argument("-i", "--input_swifter_param", required=True, help="Input Swifter parameter file to convert")
   ap.add_argument("-o", "--output_swiftest_param", required=True, help="Converted Swiftest parameter file")
   args = vars(ap.parse_args())
   inparam = args['input_swifter_param']
   outparam = args['output_swiftest_param']
   print(f"Swifter parameter is {inparam}")
   print(f"Swiftest parameter file is {outparam}")
   swio.swifter2swiftest(inparam,outparam)
