submodule (drift) s_drift_one
contains

    module procedure drift_one
    !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
    !!
    !! Perform Danby drift for one body, redoing drift with smaller substeps if original accuracy is insufficient
    !! The code has been vectorized as an elemental procedure
    !!
    !! Adapted from David E. Kaufmann's Swifter routine routine drift_one.f90
    !! Adapted from Hal Levison and Martin Duncan's Swift routine drift_one.f
     use swiftest_globals
    integer(I4B) :: i
    real(DP)     :: dttmp
    real(DP),dimension(NDIM) :: x,v
    
    x = (/posx, posy, posz/)
    v = (/vx, vy, vz/)

    call drift_dan(mu, x(:), v(:), dt, iflag)
    if (iflag /= 0) then
         dttmp = 0.1_DP * dt
         do i = 1, 10
              call drift_dan(mu, x(:), v(:), dttmp, iflag)
              if (iflag /= 0) return
         end do
    end if

    return
    end procedure drift_one

end submodule s_drift_one