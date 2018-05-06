MODULE qual_assess_contours
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! quality assess contours !!
!!
!! inputs 
!! THETA_S ::  arcmin :: Cluster Size (injected)
!! Y_inj   :: arcmin^2 CY5R500 (injected) -- all Ys are CY5R500 
!! image_contours :: 2D contours Ts & Cy5r500 :: image_contours( Ts , Cy5R00)
!! min_rs, max_rs, min_ytot, max_ytot :: beginning & end points in each dimension
!!
!! qa results array returned !!
!!
!! 1 :: HPD_2D     :: dimensionless :: HPD required to enclose the true cluster parameters
!! 2 :: HPD_Y      :: dimensionless :: HPD required to enclose the true cluster Y value (marginalise over theta)
!! 3 :: HPD_THETA  :: dimensionless :: HPD required to enclose the true cluster ThetaS value (marginalise over Y)
!! 4 :: HPD_YTHETA :: dimensionless :: HPD required to enclose the true cluster Y value (on closest grid to true theta)
!! 5 :: Ypeak     :: arcmin^2 :: Peak in Y (peak from maginalised distribution)
!! 6 :: Y_68_low  :: arcmin^2 :: lower limit of the width of HPD which encloses 68 % of the probability when theta is maginalised over
!! 7 :: Y_68_up   :: arcmin^2 :: upper limit of the width of HPD which encloses 68 % of the probability when theta is maginalised over
!! 8 :: Y_95_low  :: arcmin^2 lower limit of the width of HPD which encloses 95 % of the probability when theta is maginalised over
!! 9 :: Y_95_up   :: arcmin^2 upper limit of the width of HPD which encloses 95 % of the probability when theta is maginalised over
!!10 :: ThetaS_peak   :: arcmin :: Peak in ThetaS (peak from maginalised distribution)
!!11 :: ThetaS_68_low :: arcmin :: lower limit of the width of HPD which encloses 68 % of the probability when Y is maginalised over
!!12 :: ThetaS_68_up  :: arcmin  :: upper limit of the width of HPD which encloses 68 % of the probability when Y is maginalised over
!!13 :: ThetaS_95_low :: arcmin :: lower limit of the width of HPD which encloses 95 % of the probability when Y is maginalised over
!!14 :: ThetaS_95_up  :: arcmin :: upper limit of the width of HPD which encloses 95 % of the probability when Y is maginalised over
!!15 :: YTheta_peak   :: arcmin^2 :: Peak in Y at injected value of ThetaS (closest value on grid)
!!16 :: YTheta_68_low :: arcmin^2 :: lower limit of the width of HPD which encloses 68 % of the probability
!!17 :: YTheta_68_up  :: arcmin^2 :: upper limit of the width of HPD which encloses 68 % of the probability
!!18 :: YTheta_95_low :: arcmin^2 :: lower limit of the width of HPD which encloses 95 % of the probability
!!19 :: YTheta_95_up  :: arcmin^2 ::upper limit of the width of HPD which encloses 95 % of the probability
!!
!!20 :: Y_low_axis_prob  :: sum of all points with probabilities >= maximum probability point on the lower limit Y-axis 
!!21 :: Y_high_axis_prob :: sum of all points with probabilities >= maximum probability point on the upper limit Y-axis 
!!22 :: ThetaS_low_axis_prob  :: sum of all points with probabilities >= maximum probability point on the lower limit ThetaS-axis
!!23 :: ThetaS_high_axis_prob :: sum of all points with probabilities >= maximum probability point on the upper limit ThetaS-axis
!!
!! NB HPD set to +10 if input value is off-grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE piolibkind
  !!
  IMPLICIT NONE
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PRIVATE :: marginalise, &
       probs_from_marginal, & 
       renormalise_slice, &
       hpd_1d, &
       hpd_2d

  PUBLIC :: quality_assess_contours

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE quality_assess_contours(nscales, input_theta_s, input_cy5r500, min_rs, max_rs, min_ytot, max_ytot,&
       image_contours, qa_results) 
    INTEGER, PARAMETER :: qa_ncols=23
    !!in
    INTEGER, INTENT(in) :: nscales
    REAL(piodouble), INTENT(in) :: input_theta_s, input_cy5R500   
    REAL(piodouble), INTENT(in) :: min_rs, max_rs, min_ytot, max_ytot
    REAL(piodouble), DIMENSION(nscales,nscales), INTENT(in) :: image_contours !! image_contours( Ts , Cy5R00)
    !!out
    REAL(piodouble), DIMENSION(qa_ncols), INTENT(out) :: qa_results
    !!internal
    INTEGER :: i,j
    REAL(piodouble) :: step_yval, step_rval
    REAL(piodouble) :: hpd, peak, low68, up68, low95, up95
    REAL(piodouble), DIMENSION(4) :: edge_hpds
    REAL(piodouble), DIMENSION(nscales) :: theta_grid, y_grid 
    REAL(piodouble), DIMENSION(nscales) :: marg_theta, marg_y, slice

    !! step size in Y
    step_yval=(max_ytot-min_ytot)/REAL(nscales-1,kind=piodouble)
    !! step size in R
    step_rval=(max_rs-min_rs)/REAL(nscales-1,kind=piodouble)
    !!
    DO i=1, nscales
       theta_grid(i)=min_rs+(i-1)*step_rval !! units of arcmin
       y_grid(i)=min_ytot+(i-1)*step_yval
    ENDDO
    !!
    !! check whether input cy5r500 value lies in the range of the Y-grid AND  whether input_theta_s value lies in the range of the theta_grid
    IF((input_cy5r500.ge.y_grid(1)).and.(input_cy5r500.le.y_grid(nscales)).and.&
         (input_theta_s.ge.theta_grid(1)).and.(input_theta_s.le.theta_grid(nscales))) THEN
       !! get 2d hpd 
       CALL hpd_2d(nscales, input_theta_s, input_cy5R500, theta_grid, y_grid, image_contours, hpd)
       qa_results(1)=hpd !! HPD_2D
    ELSE
       qa_results(1)=(+10)
    ENDIF
    !!
    !! marginalise over theta => Y marginal
    CALL marginalise(0,nscales,image_contours, marg_theta)
    !! check whether input cy5r500 value lies in the range of the Y-grid
    IF((input_cy5r500.ge.y_grid(1)).and.(input_cy5r500.le.y_grid(nscales))) THEN
       CALL hpd_1d(nscales, input_cy5r500, marg_theta, y_grid, hpd)   
       qa_results(2)=hpd !! HPD_Y
    ELSE
       qa_results(2)=(+10)
    ENDIF
    !!
    !! find Y values
    CALL probs_from_marginal(nscales, marg_theta, y_grid, peak, low68, up68, low95, up95)
    !! store
    qa_results(5)=peak  !! Ypeak 
    qa_results(6)=low68 !! Y_68_low
    qa_results(7)=up68  !! Y_68_up
    qa_results(8)=low95 !! Y_95_low
    qa_results(9)=up95  !! Y_95_up
    !!
    !! marginalise over Y => Theta marginal
    CALL marginalise(1,nscales,image_contours, marg_y)
    !! check whether input_theta_s value lies in the range of the theta_grid
    IF((input_theta_s.ge.theta_grid(1)).and.(input_theta_s.le.theta_grid(nscales))) THEN
       CALL hpd_1d(nscales, input_theta_s, marg_y, theta_grid, hpd)   
       qa_results(3)=hpd !! HPD_Theta
    ELSE
       qa_results(3)=(+10)
    ENDIF
    !! find theta values
    CALL probs_from_marginal(nscales, marg_y, theta_grid, peak, low68, up68, low95, up95)
    !! store
    qa_results(10)=peak   !! ThetaS_peak
    qa_results(11)=low68 !! ThetaS_68_low
    qa_results(12)=up68  !! ThetaS_68_up
    qa_results(13)=low95 !! ThetaS_95_low
    qa_results(14)=up95  !! ThetaS_95_up
    !!
    !! get slice through input_theta_s
    CALL renormalise_slice(0, nscales, input_theta_s, theta_grid, image_contours, slice)
    CALL hpd_1d(nscales, input_cy5r500, slice, y_grid, hpd)   
    qa_results(4)=hpd !! HPD_YTheta
    !!
    !! find Y values
    CALL probs_from_marginal(nscales, slice, y_grid, peak, low68, up68, low95, up95)
    !! store
    qa_results(15)=peak  !! YThetaS_peak
    qa_results(16)=low68 !! YThetaS_68_low
    qa_results(17)=up68  !! YThetaS_68_up
    qa_results(18)=low95 !! YThetaS_95_low
    qa_results(19)=up95  !! YThetaS_95_up
    !!
    !!   
    CALL hpd_value_on_edges(nscales, image_contours, edge_hpds)
    !!store
    qa_results(20)=edge_hpds(1) !! Y_low_axis_prob 
    qa_results(21)=edge_hpds(2) !! Y_high_axis_prob 
    qa_results(22)=edge_hpds(3) !! ThetaS_low_axis_prob 
    qa_results(23)=edge_hpds(4) !! ThetaS_high_axis_prob
    !!
  END SUBROUTINE quality_assess_contours
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE hpd_2d(nscales, input_theta, input_y, theta_grid, y_grid, image_contours, hpd)
    INTEGER, INTENT(in) :: nscales
    REAL(piodouble) :: input_theta, input_y
    REAL(piodouble), DIMENSION(nscales) :: theta_grid, y_grid
    REAL(piodouble), DIMENSION(nscales,nscales), INTENT(in) :: image_contours
    REAL(piodouble), INTENT(out) :: hpd
    !!
    LOGICAL :: contains_input
    INTEGER :: i,j
    INTEGER :: peak_theta_idx, peak_y_idx
    INTEGER :: theta_idx, y_idx
    !!
    REAL(piodouble) :: maxval, total_prob, dif, minval, step
    REAL(piodouble), DIMENSION(nscales,nscales) :: probs
    !!
    !! assuming smoothly decreasing probs around the peak
    !!
    !! need to find grid position of input values
    minval=1e30
    DO i=1,nscales
       dif=ABS(theta_grid(i)-input_theta)
       IF(dif.LT.minval) THEN
          minval=dif
          theta_idx=i
       ENDIF
    ENDDO
    minval=1e30
    DO i=1,nscales
       dif=ABS(y_grid(i)-input_y)
       IF(dif.LT.minval) THEN
          minval=dif
          y_idx=i
       ENDIF
    ENDDO
    !!  
    maxval=(-1e30)
    DO i=1,nscales
       DO j=1,nscales
          IF(image_contours(i,j).GT.maxval) THEN
             maxval=image_contours(i,j)
             peak_theta_idx=i
             peak_y_idx=j
          ENDIF
          !! copy
          probs(i,j)=image_contours(i,j)
       ENDDO
    ENDDO
    !!
    total_prob=probs(peak_theta_idx,peak_y_idx)
    !!set probs->(-1) at this point
    probs(peak_theta_idx,peak_y_idx)=(-1.0)
    !!
    IF((theta_idx.eq.peak_theta_idx).and.(y_idx.eq.peak_y_idx)) then
       contains_input=.TRUE.
    ELSE
       contains_input=.FALSE.
    ENDIF
    !!
    !! check whether input values are off-grid
    !! theta
    IF((theta_idx.eq.1).or.(theta_idx.eq.nscales)) THEN
       step=theta_grid(2)-theta_grid(1)
       if(theta_idx.eq.1) then
          if((theta_grid(1)-step*0.5).gt.input_theta) then
             total_prob=1.0
             contains_input=.TRUE.
          endif
       endif
       !!
       if(theta_idx.eq.nscales) then
          if((theta_grid(nscales)+step*0.5).lt.input_theta) then
             total_prob=1.0
             contains_input=.TRUE.
          endif
       endif
    ENDIF
    !! y
    IF((y_idx.eq.1).or.(y_idx.eq.nscales)) THEN
       step=y_grid(2)-y_grid(1)
       if(y_idx.eq.1) then
          if((y_grid(1)-step*0.5).gt.input_y) then
             total_prob=1.0
             contains_input=.TRUE.
          endif
       endif
       !!
       if(y_idx.eq.nscales) then
          if((y_grid(nscales)+step*0.5).lt.input_y) then
             total_prob=1.0
             contains_input=.TRUE.
          endif
       endif
    ENDIF
    !!
    DO WHILE(.not.contains_input)
       maxval=(-1e30)
       DO i=1,nscales
          DO j=1,nscales
             IF(probs(i,j).GT.maxval) THEN
                maxval=probs(i,j)
                peak_theta_idx=i
                peak_y_idx=j
             ENDIF
          ENDDO
       ENDDO
       !!
       total_prob=total_prob+probs(peak_theta_idx,peak_y_idx)
       !!set probs->(-1) at this point
       probs(peak_theta_idx,peak_y_idx)=(-1.0)
       !!
       IF((theta_idx.eq.peak_theta_idx).and.(y_idx.eq.peak_y_idx)) contains_input=.TRUE.
       !!
    ENDDO
    !!
    hpd=total_prob
    !!
  END SUBROUTINE hpd_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE hpd_1d(nscales, input_value, marginal_probs, marginal_quants, hpd)
    INTEGER, INTENT(in) :: nscales
    REAL(piodouble) :: input_value
    REAL(piodouble), DIMENSION(nscales), INTENT(in) :: marginal_probs, marginal_quants
    REAL(piodouble), INTENT(out) :: hpd
        !!
    LOGICAL :: contains_input, lower_bound_hit, upper_bound_hit, trouble
    INTEGER :: i, peak_idx, low_idx, up_idx, test_low, test_up
    REAL(piodouble) :: maxval, total_prob, stepsize
    REAL(piodouble) :: lower_lim, upper_lim
    !!
    maxval=(-1e30)
    DO i=1,nscales
       IF(marginal_probs(i).GT.maxval) THEN
          maxval=marginal_probs(i)
          peak_idx=i
       ENDIF
    ENDDO
    !!
    stepsize=marginal_quants(nscales)-marginal_quants(nscales-1)
    lower_lim=marginal_quants(peak_idx)-stepsize*0.5
    upper_lim=marginal_quants(peak_idx)+stepsize*0.5
    !!
    !! assuming smoothly decreasing probs on both sides of the peak
    !!
    total_prob=marginal_probs(peak_idx)
    !!
    IF((input_value.le.upper_lim).and.(input_value.ge.lower_lim)) then
       contains_input=.TRUE.
    ELSE
       contains_input=.FALSE.
    ENDIF
    !!
    lower_bound_hit=.FALSE.
    upper_bound_hit=.FALSE.
    trouble=.FALSE.
    low_idx=peak_idx
    up_idx=peak_idx
    !! check bounds
    IF(low_idx.EQ.1) lower_bound_hit=.TRUE. 
    IF(up_idx.EQ.nscales) upper_bound_hit=.TRUE. 
    !!
    !! check if input is off-grid
    IF(lower_bound_hit) THEN
       IF((marginal_quants(1)-stepsize*0.5).gt.input_value) THEN
          total_prob=1.0
          contains_input=.TRUE.
       ENDIF
    ENDIF
    IF(upper_bound_hit) THEN
       IF((marginal_quants(nscales)+stepsize*0.5).lt.input_value) THEN
          total_prob=1.0
          contains_input=.TRUE.
       ENDIF
    ENDIF
    !!
    DO WHILE((.not.contains_input).and.(.not.trouble))
       IF((.NOT.lower_bound_hit).AND.(.NOT.upper_bound_hit)) THEN
          !! advance 
          test_low=low_idx-1
          test_up=up_idx+1
          !!
          IF(marginal_probs(test_low).GT.marginal_probs(test_up)) THEN
             !! add prob of marginal_probs(test_low)
             total_prob=total_prob+marginal_probs(test_low)
             !!advance low_idx
             low_idx=low_idx-1
             !! check bounds
             IF(low_idx.EQ.1) lower_bound_hit=.TRUE. 
             !!
          ELSE
             !! add prob of marginal_probs(test_up)
             total_prob=total_prob+marginal_probs(test_up)
             !!advance up_idx
             up_idx=up_idx+1
             !! check bounds
             IF(up_idx.EQ.nscales) upper_bound_hit=.TRUE. 
             !!
          ENDIF
       ELSE
          !! one of the boundaries has been hit
          IF(.NOT.lower_bound_hit) THEN
             !!advance low_idx
             low_idx=low_idx-1
             !! add prob
             total_prob=total_prob+marginal_probs(low_idx)
             !! check bounds
             IF(low_idx.EQ.1) lower_bound_hit=.TRUE. 
          ENDIF
          !!
          IF(.NOT.upper_bound_hit) THEN
             !! advance up_idx
             up_idx=up_idx+1
             !! add prob
              total_prob=total_prob+marginal_probs(up_idx)
             !! check bounds
             IF(up_idx.EQ.nscales) upper_bound_hit=.TRUE. 
          ENDIF
          !!
          IF(lower_bound_hit.and.upper_bound_hit) THEN
             print*,'total_prob=',total_prob,' rather than 1'
             trouble=.TRUE.
          ENDIF
       ENDIF
       !!
       !! check to see whether the input value now lies withing the area of summed probs.
       !!
       lower_lim=marginal_quants(low_idx)-stepsize*0.5
       upper_lim=marginal_quants(up_idx)+stepsize*0.5
       IF((input_value.le.upper_lim).and.(input_value.ge.lower_lim)) contains_input=.TRUE.
       !!
    ENDDO
    !!
    hpd=total_prob
    !!
  END SUBROUTINE hpd_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE renormalise_slice(ndim, nscales, slice_point, slice_grid, image_contours, slice)
    !! ndim=0 slice point on rows :: slice runs through cols
    !! ndim=1 slice point on cols :: slice runs through rows
    !! image_contours(rows,cols)
    INTEGER, INTENT(in) :: ndim, nscales
    REAL(piodouble), INTENT(in) :: slice_point
    REAL(piodouble), DIMENSION(nscales), INTENT(in) :: slice_grid
    REAL(piodouble), DIMENSION(nscales,nscales), INTENT(in) :: image_contours !! image_contours( Ts , Cy5R00)
    REAL(piodouble), DIMENSION(nscales), INTENT(out) :: slice
    !!
    INTEGER :: i,j
    INTEGER :: min_idx
    REAL(piodouble) :: norm, dif, minval
    !!
    !! find closest point in grid to slice_point
    minval=1e30
    DO i=1,nscales
       dif=ABS(slice_grid(i)-slice_point)
       IF(dif.LT.minval) THEN
          minval=dif
          min_idx=i
       ENDIF
    ENDDO
    !! get slice     
    norm=0.0
    IF(ndim.EQ.0) THEN
       DO i=1,nscales
          slice(i)=image_contours(min_idx,i)
          norm=norm+slice(i)
       ENDDO
    ENDIF
    IF(ndim.EQ.1) THEN
       DO i=1,nscales
          slice(i)=image_contours(i,min_idx)
          norm=norm+slice(i)
       ENDDO
    ENDIF
    !! renormalise
    DO i=1,nscales
       slice(i)=slice(i)/norm
    ENDDO
    !!
  END SUBROUTINE renormalise_slice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE marginalise(ndim, nscales, image_contours, marginal)
    !! ndim=0 marginalise over rows
    !! ndim=1 marginalise over cols
    !! image_contours(rows,cols)
    INTEGER, INTENT(in) :: ndim, nscales
    REAL(piodouble), DIMENSION(nscales,nscales), INTENT(in) :: image_contours !! image_contours( Ts , Cy5R00)
    REAL(piodouble), DIMENSION(nscales), INTENT(out) :: marginal
    !!
    INTEGER :: i,j
    REAL(piodouble) :: sum
    !!
    IF(ndim.EQ.0) THEN
       DO i=1,nscales
          sum=0.0
          DO j=1,nscales
             sum=sum+image_contours(j,i) !! add up rows
          ENDDO
          marginal(i)=sum
       ENDDO
    ENDIF
    !!
    IF(ndim.EQ.1) THEN
       DO i=1,nscales
          sum=0.0
          DO j=1,nscales
             sum=sum+image_contours(i,j) !! add up cols
          ENDDO
          marginal(i)=sum
       ENDDO
    ENDIF
    !!
  END SUBROUTINE marginalise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE probs_from_marginal(nscales, marginal_probs, marginal_quants,&
        peak, low68, up68, low95, up95)
    REAL(piodouble), PARAMETER :: first_limit=0.68
    REAL(piodouble), PARAMETER :: second_limit=0.95
    !!
    INTEGER, INTENT(in) :: nscales
    REAL(piodouble), DIMENSION(nscales), INTENT(in) :: marginal_probs, marginal_quants
    REAL(piodouble), INTENT(out) :: peak, low68, up68, low95, up95
    !!
    LOGICAL :: first_limit_found, lower_bound_hit, upper_bound_hit, trouble
    INTEGER :: i, peak_idx, low_idx, up_idx, test_low, test_up
    REAL(piodouble) :: maxval, total_prob
    !!
    maxval=(-1e30)
    DO i=1,nscales
       IF(marginal_probs(i).GT.maxval) THEN
          maxval=marginal_probs(i)
          peak_idx=i
       ENDIF
    ENDDO
    !!
    !! have value at peak
    peak=marginal_quants(peak_idx)
    !!
    !! assuming smoothly decreasing probs on both sides of the peak
    !!
    total_prob=marginal_probs(peak_idx)
    first_limit_found=.FALSE.
    lower_bound_hit=.FALSE.
    upper_bound_hit=.FALSE.
    trouble=.FALSE.
    low_idx=peak_idx
    up_idx=peak_idx
    !! check bounds
    IF(low_idx.EQ.1) lower_bound_hit=.TRUE. 
    IF(up_idx.EQ.nscales) upper_bound_hit=.TRUE. 
    DO WHILE((total_prob.LT.second_limit).and.(.not.trouble))
       IF((.NOT.lower_bound_hit).AND.(.NOT.upper_bound_hit)) THEN
          !! advance 
          test_low=low_idx-1
          test_up=up_idx+1
          !!
          IF(marginal_probs(test_low).GT.marginal_probs(test_up)) THEN
             !! add prob of marginal_probs(test_low)
             total_prob=total_prob+marginal_probs(test_low)
             !!advance low_idx
             low_idx=low_idx-1
             !! check bounds
             IF(low_idx.EQ.1) lower_bound_hit=.TRUE. 
             !!
          ELSE
             !! add prob of marginal_probs(test_up)
             total_prob=total_prob+marginal_probs(test_up)
             !!advance up_idx
             up_idx=up_idx+1
             !! check bounds
             IF(up_idx.EQ.nscales) upper_bound_hit=.TRUE. 
             !!
          ENDIF
       ELSE
          !! one of the boundaries has been hit
          IF(.NOT.lower_bound_hit) THEN
             !!advance low_idx
             low_idx=low_idx-1
             !! add prob
             total_prob=total_prob+marginal_probs(low_idx)
             !! check bounds
             IF(low_idx.EQ.1) lower_bound_hit=.TRUE. 
          ENDIF
          !!
          IF(.NOT.upper_bound_hit) THEN
             !! advance up_idx
             up_idx=up_idx+1
             !! add prob
              total_prob=total_prob+marginal_probs(up_idx)
             !! check bounds
             IF(up_idx.EQ.nscales) upper_bound_hit=.TRUE. 
          ENDIF
          !!
          IF(lower_bound_hit.and.upper_bound_hit) THEN
             print*,'total_prob=',total_prob,' rather than 1'
             trouble=.TRUE.
          ENDIF
       ENDIF
       !!
       IF(.NOT.first_limit_found) THEN
          IF(total_prob.GE.first_limit) THEN
             low68=marginal_quants(low_idx)
             up68=marginal_quants(up_idx)
             first_limit_found=.TRUE.
          ENDIF
       ENDIF
    ENDDO
    !!
    low95=marginal_quants(low_idx)
    up95=marginal_quants(up_idx)
    !!
  END SUBROUTINE probs_from_marginal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE  hpd_value_on_edges(nscales, image_contours, edge_hpds)
    !!
    INTEGER, INTENT(in) :: nscales
    REAL(piodouble), DIMENSION(nscales,nscales), INTENT(in) :: image_contours !! image_contours( Ts , Cy5R00)
    REAL(piodouble), DIMENSION(4), INTENT(out) :: edge_hpds
    !!
    INTEGER :: i,j, found_count
    INTEGER :: peak_theta_idx, peak_y_idx
    INTEGER :: lowerYaxis_idx, upperYaxis_idx
    INTEGER :: lowerTs_axis_idx, upperTs_axis_idx
    REAL(piodouble) :: maxval, total_prob
    REAL(piodouble), DIMENSION(4) :: edge_max_probs
    REAL(piodouble), DIMENSION(nscales,nscales) :: probs
    !!
    !! initialise
    found_count=0
    total_prob=0.0
    edge_max_probs=0.0
    DO i=1,nscales
       DO j=1,nscales
          probs(i,j)=image_contours(i,j)
       ENDDO
    ENDDO
    !!
    lowerYaxis_idx=1
    upperYaxis_idx=1
    lowerTs_axis_idx=nscales
    upperTs_axis_idx=nscales
    !!
    DO i=1,nscales
       !! lower Yval axis :: index=1
       if(edge_max_probs(1).lt.image_contours(i,1)) then
          edge_max_probs(1)=image_contours(i,1) !! thetaS varies / Y5r500 const
          lowerYaxis_idx=i
       endif
       !! upper Yval axis :: index=2
       if(edge_max_probs(2).lt.image_contours(i,nscales)) then
          edge_max_probs(2)=image_contours(i,nscales) 
          upperYaxis_idx=i
       endif
       !! lower ThetaS axis :: index=3
       if(edge_max_probs(3).lt.image_contours(1,i)) then
          edge_max_probs(3)=image_contours(1,i) !! Y5r500 varies / thetas Const
          lowerTs_axis_idx=i
       endif
       !! upper ThetaS axis :: index=4 
       if(edge_max_probs(4).lt.image_contours(nscales,i)) then
          edge_max_probs(4)=image_contours(nscales,i) 
          upperTs_axis_idx=i
       endif
    ENDDO
    !! 
    !! assuming smoothly decreasing probs around one peak
    !!
    DO WHILE(found_count.lt.4)
       maxval=(-1e30)
       DO i=1,nscales
          DO j=1,nscales
             IF(probs(i,j).GT.maxval) THEN
                maxval=probs(i,j)
                peak_theta_idx=i
                peak_y_idx=j
             ENDIF
          ENDDO
       ENDDO
       !!
       total_prob=total_prob+probs(peak_theta_idx,peak_y_idx)
       !!set probs->(-1) at this point
       probs(peak_theta_idx,peak_y_idx)=(-1.0)
       !!
       if((peak_y_idx.eq.1).and.(peak_theta_idx.eq.lowerYaxis_idx)) then
          !! lower Yval axis :: index=1
          edge_hpds(1)=total_prob
          found_count=found_count+1
       endif
       if((peak_y_idx.eq.nscales).and.(peak_theta_idx.eq.upperYaxis_idx)) then
          !! upper Yval axis :: index=2
          edge_hpds(2)=total_prob
          found_count=found_count+1
       endif
       if((peak_theta_idx.eq.1).and.(peak_y_idx.eq.lowerTs_axis_idx)) then
          !! lower ThetaS axis :: index=3
          edge_hpds(3)=total_prob
          found_count=found_count+1
       endif
       if((peak_theta_idx.eq.nscales).and.(peak_y_idx.eq.upperTs_axis_idx)) then
          !! upper ThetaS axis :: index=4 
          edge_hpds(4)=total_prob
          found_count=found_count+1
       endif
    ENDDO
    !!
  END SUBROUTINE hpd_value_on_edges
  !!
END MODULE qual_assess_contours
