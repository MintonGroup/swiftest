submodule (util) s_util_version
contains
   module procedure util_version
   !! author: David A. Minton
   !!
   !! Print program version information to terminale
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: util_version.f90
   use swiftest

   write(*, 200) VERSION_NUMBER
200 format(/, "************* Swiftest: Version ", f3.1, " *************", //, &
         "Based off of Swifter:", //,                                         &
         "Authors:", //,                                                      &
         "    The Purdue University Swiftest Development team ", /,           &
         "    Lead by David A. Minton ", /,                                   &
         "    Single loop blocking by Jacob R. Elliott", /,                   &
         "    Fragmentation by Carlisle A. Wishard and", //,                  &
         "    Jennifer L. L. Poutplin                 ", //,                  &
         "Please address comments and questions to:", //,                     &
         "    David A. Minton", /,                                            &
         "    Department Earth, Atmospheric, & Planetary Sciences ",/,        &
         "    Purdue University", /,                                          &
         "    550 Stadium Mall Drive", /,                                     &
         "    West Lafayette, Indiana 47907", /,                              &
         "    765-250-8034 ", /,                                              &
         "    daminton@purdue.edu", /,                                        &
         "Special thanks to Hal Levison and Martin Duncan for the original",/,&
         "SWIFTER and SWIFT codes that made this possible.", //,              &
         "************************************************", /)


100 FORMAT(/,  "************* SWIFTER: Version ", F3.1, " *************", //, &
               "Authors:", //,                                                &
               "    Martin Duncan: Queen's University", /,                    &
               "    Hal Levison  : Southwest Research Institute", //,         &
               "Please address comments and questions to:", //,               &
               "    Hal Levison or David Kaufmann", /,                        &
               "    Department of Space Studies", /,                          &
               "    Southwest Research Institute", /,                         &
               "    1050 Walnut Street, Suite 400", /,                        &
               "    Boulder, Colorado  80302", /,                             &
               "    303-546-0290 (HFL), 720-240-0119 (DEK)", /,               &
               "    303-546-9687 (fax)", /,                                   &
               "    hal@gort.boulder.swri.edu (HFL)", /,                      &
               "    kaufmann@boulder.swri.edu (DEK)", //,                     &
               "************************************************", /)

     return

   end procedure util_version
end submodule s_util_version
