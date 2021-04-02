#!/bin/bash
for file_out in */*.f90; do
   file_in="../../swifter-omp/$file_out";
   desc=$(grep "Description" $file_in  | sed "s/!  Description : //")
   sed -i "" "s/Compute Hill sphere radii of massive bodie/$desc/" $file_out
done

