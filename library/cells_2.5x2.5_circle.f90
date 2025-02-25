! Processing data from electrochemical workstation to obtain 
! cell properties
! by zito, 2018/12/17
!
! some bugs fixed
! Parameter 'n_slop_average' is added
! by zito, 2018/12/21
!
! slope of Rs are calculated by three points
! by zito, 2018/12/24 

      PROGRAM cell_property
      IMPLICIT NONE
      INTEGER, PARAMETER       :: sr=4,sc=8,dr=8,dc=16
! n_slop_average*2 + 1 is the number of data used to calculate slope of
! linear regression line
      INTEGER, PARAMETER       :: n_slop_average=50
      REAL(KIND=dr), PARAMETER :: area=0.19625d0
      CHARACTER(LEN=40)        :: filename
      CHARACTER(LEN=79)        :: buffer
      CHARACTER(LEN=9)         :: read_flag
      INTEGER                  :: i_read,i_file,i,n_row
      REAL(KIND=dr)            :: voltage,current,current_min, &
                                  v_current_min,v_zero,current_v_zero, &
                                  v_x_c,pce,ff,temp,slope_rsh,slope_rs,&
                                  voltage_ave, current_ave, sum_xy, &
                                  sum_x2
      REAL(KIND=dr),DIMENSION(-n_slop_average:n_slop_average) :: & 
                                  voltage_col,current_col
!
      filename="nonefile"
      OPEN(UNIT=7,FILE="list.txt",STATUS="OLD")
! Table title
        WRITE(*,'(A)') "File_number     Voc     Jsc     PCE     FF      &
          Rsh   Rs"
      DO WHILE(.TRUE.)
        READ(UNIT=7,FMT="(A)") filename
        IF (TRIM(ADJUSTL(filename)) .EQ. "EOF") EXIT
!      
        OPEN(UNIT=77,FILE=filename,STATUS="OLD")
!        WRITE(*,'(2A)') "The file name is ", filename
! discard useless text
        read_flag="reading  "
        DO WHILE (read_flag .NE. "Potential")
          READ(77,"(A9)") read_flag
        ENDDO
        READ(77,"(A9)") read_flag
! counting total number of rows
        i_read=0
        n_row=0
        DO WHILE (i_read .EQ. 0)
          n_row=n_row+1
          READ(77, FMT="(A79)", IOSTAT=i_read) buffer
!          WRITE(*,*) buffer, n_row
        ENDDO
!
        REWIND(77)
        read_flag="reading  "
        DO WHILE (read_flag .NE. "Potential")
          READ(77,"(A9)") read_flag
!          WRITE(*,*) read_flag
        ENDDO
        READ(77,"(A9)") read_flag
! read voltage and current, and find the most approching zero current
! value and the corresponding voltage vaule
        current_min=100.d0
        pce=1.d-20
        DO i=1,n_row-1
          READ(77,*) voltage,current
!          WRITE(*,*) voltage,current
! find Voc
          IF (ABS(current) .LE. ABS(current_min)) THEN
            current_min=current
            v_current_min=voltage
          ENDIF
! find Jsc
          IF (voltage .EQ. 0.d0) THEN
            v_zero=voltage
            current_v_zero=current
          ENDIF
! PCE
          v_x_c=-voltage*current
!          WRITE(*,*) voltage, current, v_x_c
          IF (v_x_c .GT. pce ) THEN
            pce=v_x_c
          ENDIF
        ENDDO
        pce=pce*10/area
! FF
        ff=ABS(pce/v_current_min/(current_v_zero*1000/area))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!                 Calculate Rsh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! turn back to the begining of the data
        REWIND(77)
        read_flag="reading  "
        DO WHILE (read_flag .NE. "Potential")
          READ(77,"(A9)") read_flag
!          WRITE(*,*) read_flag
        ENDDO
        READ(77,"(A9)") read_flag
! now read data
        voltage=1.d0
        current=1.d0
        DO WHILE (voltage .NE. 0.d0) 
          READ(77,*) voltage,current
        ENDDO
!        WRITE(*,*) voltage,current
! backspace n_slop_average+1 lines
        DO i=1,n_slop_average+1
          BACKSPACE(77)
        ENDDO
! average n_slop_average*2+1 data
        voltage_col(:)=0.d0
        current_col(:)=0.d0
        voltage_ave=0.d0
        current_ave=0.d0
        sum_xy=0.d0
        sum_x2=0.d0
        slope_rsh=0.d0
        DO i = -n_slop_average, n_slop_average, 1
          READ(77,*) voltage_col(i),current_col(i)
          voltage_ave=voltage_ave+voltage_col(i)
          current_ave=current_ave+current_col(i)
!          WRITE(*,*) voltage_col(i),current_col(i)
        ENDDO
        voltage_ave=voltage_ave/100
        current_ave=current_ave/100
        DO i = -n_slop_average, n_slop_average, 1
          sum_xy = sum_xy + (voltage_col(i)-voltage_ave) * &
                   (current_col(i)-current_ave)
          sum_x2 = sum_x2 + (voltage_col(i)-voltage_ave) ** 2
        ENDDO
        slope_rsh= sum_x2 / sum_xy
!        WRITE(*,'(A11,F14.5)') "RSH     :  ", slope_rsh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 Calculate Rs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! turn back to the begining of the data
        REWIND(77)
        read_flag="reading  "
        DO WHILE (read_flag .NE. "Potential")
          READ(77,"(A9)") read_flag
!          WRITE(*,*) read_flag
        ENDDO
        READ(77,"(A9)") read_flag
! now read data
        voltage=1.d0
        current=1.d3
        DO WHILE (current .NE. current_min)
          READ(77,*) voltage,current
        ENDDO
!        WRITE(*,*) voltage,current 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! backspace n_slop_average+1 lines
!        DO i=1,n_slop_average-45+1
!          BACKSPACE(77)
!        ENDDO
! average n_slop_average*2+1 data
!       voltage_col(:)=0.d0
!       current_col(:)=0.d0
!       voltage_ave=0.d0
!       current_ave=0.d0
!       sum_xy=0.d0
!       sum_x2=0.d0
!       slope_rs=0.d0
!       DO i = -n_slop_average+45, n_slop_average-45, 1
!         READ(77,*) voltage_col(i),current_col(i)
!         voltage_ave=voltage_ave+voltage_col(i)
!         current_ave=current_ave+current_col(i)
!          WRITE(*,*) voltage_col(i),current_col(i)
!       ENDDO
!       voltage_ave=voltage_ave/100
!       current_ave=current_ave/100
!       DO i = -n_slop_average+45, n_slop_average-45, 1
!         sum_xy = sum_xy + (voltage_col(i)-voltage_ave) * &
!                  (current_col(i)-current_ave)
!         sum_x2 = sum_x2 + (voltage_col(i)-voltage_ave) ** 2
!       ENDDO
!       slope_rs= sum_x2 / sum_xy * area
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! backspace  one line
        BACKSPACE(77)
        voltage_col(:)=0.d0
        current_col(:)=0.d0
        slope_rs=0.d0
        DO i = -1, 1, 1
          READ(77,*) voltage_col(i),current_col(i)
        ENDDO
        slope_rs = (current_col(1)-current_col(-1)) / ((voltage_col(1)) &
                   - voltage_col(-1))
        slope_rs = slope_rs * area 
! data output
        filename=TRIM(ADJUSTL(filename))
        WRITE(*,'(A10,6F16.5)') filename, ABS(v_current_min), &
                           ABS(current_v_zero*1000/area),     &
                           ABS(pce), ABS(ff*100), ABS(slope_rsh),    &
                           ABS(slope_rs)
! close data file
        CLOSE(77)

      ENDDO
      CLOSE(7)
!
      END
