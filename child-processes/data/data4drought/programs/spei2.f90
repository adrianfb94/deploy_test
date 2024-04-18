!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c       Standardized Precipitation Evapotranspiration Index.
!        This program is based on SPEI c++ algorit as well as
!        the SPI program. This was adapted to read and write
!        GrADS files, being more ease the work amny people need
!        to develop
!
!
!c       spei_grads should be executed after thorthwaite and waterbalance
!        programs, so it should be not possible to use independently.
!
!        A matrix of 100X100 elements will allowed for 250 times
!        
!        This programe compute the SPEI for 9 different periods. This means
!        1, 3, 6 12, 18, 24 36 and 48 months
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
program speiprg
parameter (inlen=7)
character(len=:),allocatable::ifile,ifile2,ofile,result_file 
character(len=:),allocatable::arg
real,allocatable,dimension(:,:,:)::tarray,spiout,wtbalance,idx,droughts
real,allocatable,dimension(:,:)::totumbral,conseumbral,longr,rndrought,prcumbral
real,allocatable,dimension(:)::idx2
real rumbral,sumbral
real miss,umbral,lat,dlat
real dlim(7)
real logLogisticParams,standardGaussianInvCDF,around,rtlen2,rspilen
integer line(20),ix,iy,nyears,yr,inx,iny,ncd,tlen,tlen2,spilen,totlen,eta,etr,iarg,ndrou

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Lectura de los parametros pasados en linea de comando

!ichars=500! esto era para pasar a allocate de los nombre de files pero no funciono
allocate(character(len=500)::arg)
allocate(character (len=500) ::ifile,ifile2,ofile,result_file)! esto es para leer lasgos nombres de files
n=iargc() !n almacena en numero total de parametros en la cmdline

 call getarg(1,arg)
 ifile=arg !nombre fichero con datos de periodo de referencia
 call getarg(2,arg)
 ifile2=arg !nombre fichero de datos de lluvia. Solo valido para el SPEI
 call getarg(3,arg)
 ofile=arg            !nombre de fichero salida del spi
 call getarg(4,arg)
 result_file=arg      !nombre fichero con datos de conteo de sequias acorde al umbral
do iarg=5,n-9
 call getarg(iarg,arg)
 read(arg,*)line(iarg-4) !parametros tlen,nmon,inx,iny,miss,tlen2,lenspi,ndrou
enddo
 call getarg(15,arg)
 read(arg,*) lat     !lat es la latitud inicial del fichero de datos de entrada
 call getarg(16,arg)
 read(arg,*) dlat    !dlat es el incremento de la latitud en el fichero
 call getarg(17,arg)
 read(arg,*)line(11)
 call getarg(18,arg)
 read(arg,*)line(12)
 call getarg(19,arg)
 read(arg,*)line(13)
 call getarg(n,arg)
 read(arg,*) umbral     !este es el unico parametro que se lee fuera del ciclo por ser real
 call getarg(20,arg)
 read(arg,*) umbral     !este es el unico parametro que se lee fuera del ciclo por ser real
 call getarg(21,arg)
 read(arg,*) dlim(1)     
  call getarg(22,arg)
 read(arg,*) dlim(2)     
 call getarg(23,arg)
 read(arg,*) dlim(3)    
 dlim(4)=abs(dlim(3))
                        !si se quiere annadir mas parametros, se recomienda sea antes de este


! Definition of variables from cmdline parameters
totlen=line(1) ! periodo total (meses*annos) del fichero de datos
tlen=line(2)   ! periodo de la linea base (meses*annos)
itr=line(3)    ! indice de inicio de periodo linea base
etr=line(4)    ! indice de final de periodo linea base
tlen2=line(5)  ! periodo de analisis (meses*annos)
ita=line(6)    ! indice de inicio de periodo analisis
eta=line(7)    ! indice de final de periodo analisis
nmon=line(8)   ! numero de meses [ahora es 12] NO TIENE MUCHA UTILIDAD
inx=line(9)    ! numero de puntos por las x (longitud)
iny=line(10)   ! numero de puntos por las y (latitud)
miss=line(11)*1. !valor missing
spilen=line(12)  !numero de meses para calcular el spi
ndrou=line(13) ! cantidad de meses limite para contar periodos validos de sequia
nyears=tlen/nmon
nyears2=tlen2/nmon
IRL=inx*iny*4   ! numero de bytes para el registro
! end of variable definition from cmdline parameters

rtlen2=tlen2*1.
rspilen=spilen*1.

allocate(tarray(inx,iny,totlen)) ! allocate the array for temperature

!definir un cmd line para leer evapotrans directamente y no calcular ETP. Entonces, aqui hay que poner un if
! y leer el file de ETP
!***********************************************************************************************************

open(31,file=ifile,FORM='UNFORMATTED',access='direct',RECL=IRL) !temperature file
IREC=1
do i=1,totlen
  read(31,REC=IREC) tarray(:,:,i)!read all the records (tarray(:,:,totlen)) as calculated by GrADS script
  IREC=IREC+1
enddo

write(*,*)'tarray',minval(tarray), maxval(tarray)


ifirst=1
last=tlen2 !here tlen2 is the lenght of the assessment record
close(31)


allocate(wtbalance(inx,iny,totlen))

do iy=1,iny
 do ix=1,inx
  do i=1,totlen
   if(tarray(ix,iy,i).ne.miss)then
    wtbalance(ix,iy,i)=around(tarray(ix,iy,i),1)
   else
    wtbalance(ix,iy,i)=miss
   endif
  enddo
 enddo
enddo

write(*,*)'wtbalance', minval(wtbalance), maxval(wtbalance)
deallocate(tarray)

allocate(spiout(inx,iny,last)) ! allocate spiout variable and other drought accountability variables

allocate(idx(inx,iny,totlen),idx2(tlen2))

allocate(longr(inx,iny),rndrought(inx,iny),totumbral(inx,iny),prcumbral(inx,iny))

allocate(droughts(inx,iny,inlen))

do iy=1,iny
 do ix=1,inx
  if(wtbalance(ix,iy,itr).ne.miss.and.wtbalance(ix,iy,ita).ne.miss)then
   call speilog(wtbalance(ix,iy,itr:etr),wtbalance(ix,iy,ita:eta),tlen,spilen,spiout(ix,iy,ifirst:last),miss,tlen2)
  !  write(*,*)'spiout(ix,iy,ifirst)', ix,iy,spiout(ix,iy,ifirst)

   !   counting some total spi values below umbral
    totumbral(ix,iy)=0 !total of spi below umbral  
    longr(ix,iy)=0 ! maximum consecutive spi values below umbral
    rumbral=0 ! temporal counting variable for maximum consecutive spi values
    rndrought(ix,iy)=0
    sumbral=0
    do i=ifirst,last
       if(spiout(ix,iy,i).le.umbral.and.spiout(ix,iy,i).ne.miss)then
        totumbral(ix,iy)=totumbral(ix,iy)+1
        rumbral=rumbral+1
        sumbral=sumbral+1
       else
        if(sumbral.ge.ndrou) then
          rndrought(ix,iy)=rndrought(ix,iy)+1 ! number of drought events greater than ndrou param
          sumbral=0
        endif
        if(rumbral.gt.longr(ix,iy))then ! determining the maximum consecutive spi period
         longr(ix,iy)=rumbral
         rumbral=0
        endif
       endif
       !this section we count the droguth categories.
       !Four categories are defined by the user throughout the parameter list
       !these are dlim(1:4) for extreme, severe, moderate and normal
       !the wet categories have fixed values in this version
       !but can be changed in the future
       if(spiout(ix,iy,i).le.dlim(1))droughts(ix,iy,1)=droughts(ix,iy,1)+1
       if(spiout(ix,iy,i).gt.dlim(1).and.spiout(ix,iy,i).le.dlim(2))droughts(ix,iy,2)=droughts(ix,iy,2)+1
       if(spiout(ix,iy,i).gt.dlim(2).and.spiout(ix,iy,i).le.dlim(3))droughts(ix,iy,3)=droughts(ix,iy,3)+1
       if(spiout(ix,iy,i).gt.dlim(3).and.spiout(ix,iy,i).le.dlim(4))droughts(ix,iy,4)=droughts(ix,iy,4)+1
       if(spiout(ix,iy,i).ge.1.00.and.spiout(ix,iy,i).le.1.49)droughts(ix,iy,5)=droughts(ix,iy,5)+1
       if(spiout(ix,iy,i).ge.1.50.and.spiout(ix,iy,i).le.1.99)droughts(ix,iy,6)=droughts(ix,iy,6)+1
       if(spiout(ix,iy,i).ge.2.0)droughts(ix,iy,7)=droughts(ix,iy,7)+1
    enddo
    prcumbral(ix,iy)=(around(totumbral(ix,iy),1)/(rtlen2-rspilen))*100
    do id=1,7
     droughts(ix,iy,id)=(around(droughts(ix,iy,id),1)/(rtlen2-rspilen))*100
    enddo
   else
     spiout(ix,iy,ifirst:last)=miss
     totumbral(ix,iy)=miss
     rndrought(ix,iy)=miss
     droughts(ix,iy,1:7)=miss
     prcumbral(ix,iy)=miss
     longr(ix,iy)=miss
   endif
 enddo
enddo
deallocate(wtbalance)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The speiout matrix is then saved into speiout.grd file wich will
!   be call by GrADS directally after the program stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c  Writing SPI values (spiout variable)

write(*,*)minval(spiout), maxval(spiout)


open(49,file=ofile,access='direct',recl=IRL)
ir=1
do i=ifirst, last
  write(49,rec=ir) spiout(:,:,i) !se escriben los valores de spi para el area y periodo definido
  ir=ir+1
enddo
close(49)
deallocate(spiout)


ofile=trim(result_file)//'_1.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) totumbral(:,:)
close(49)
deallocate(totumbral)

ofile=trim(result_file)//'_2.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) longr(:,:)
close(49)

ofile=trim(result_file)//'_3.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) rndrought(:,:)
close(49)

ofile=trim(result_file)//'_4.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) prcumbral(:,:)
close(49)

deallocate(longr,rndrought,prcumbral)

ofile=trim(result_file)//'d_1.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) droughts(:,:,1)
close(49)

ofile=trim(result_file)//'d_2.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) droughts(:,:,2)
close(49)

ofile=trim(result_file)//'d_3.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) droughts(:,:,3)
close(49)

ofile=trim(result_file)//'d_4.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) droughts(:,:,4)
close(49)

ofile=trim(result_file)//'d_5.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) droughts(:,:,5)
close(49)

ofile=trim(result_file)//'d_6.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) droughts(:,:,6)
close(49)

ofile=trim(result_file)//'d_7.grb'
open(49,file=ofile,access='direct',recl=IRL)
write(49,rec=1) droughts(:,:,7)
close(49)

deallocate(droughts)


 
 stop
 1000 format(a80)
 1001 format(i4,i3,8f9.2)
      end


subroutine speilog(idx,idx2,n,season,speiSeries,rmiss,nf)
 integer n,ik,nf,ki,a
 integer i, j, k, nSeason,nSeason2,season,fSeason
 real seasonSeries2(nf),seasonSeries(n+nf),beta(3),logLogisticParams(48,3),rmiss
 real dataSeries(n),speiSeries(nf),logpara(3),logLogisticCDF
 real dataSeries2(nf),idx(n),idx2(nf)

!  do i=1,3
!   write(*,*)logpara(i), beta(i)
! end do
! write(*,*)' '


 call cumsum(season,idx,idx2,dataSeries,dataSeries2,rmiss,n,nf)
do j=1,12
!!Tomando solo baseline para construir la funcion. 
	ki=1
	do i=season+j-1,n,12   !delete the missing values.
     if(dataSeries(i) .ne. rmiss)then
      seasonSeries(ki)=dataSeries(i)
     endif
     ki=ki+1
    enddo
	nSeason=ki-1
    call ssort(seasonSeries, nSeason)					!sort upward the data serie
!     call pwm(seasonSeries, nSeason,beta, -0.35, 0.,0.)
    call pwm(seasonSeries, nSeason,beta, 0., 0.,0.)
  !  write(*,*) nSeason,beta(1:3)
    call logLogisticFit(beta, logpara)
    
    do i=1,season-1
		speiSeries(i) = rmiss
    end do
             
    do i=season+j-1,nf,12
          
     if(dataSeries2(i).ne.rmiss)then
      ! write(*,*)'estoy en dataSeries no miss (old) i=', i,dataSeries2(i), logpara
      speiSeries(i) = logLogisticCDF(dataSeries2(i), logpara)
      ! write(*,*)'estoy en dataSeries no miss (new) i=', i,dataSeries2(i)

    else
      speiSeries(i) = rmiss
     endif
     
    enddo
   enddo

  !  do i=1,3
  !   write(*,*)logpara(i), beta(i)
  ! end do
  ! write(*,*)' '

   return
end

subroutine cumsum(nr,p,p2,idx,idx2,miss,iyrs,iyrs2)
  real p(iyrs),p2(iyrs2),idx(iyrs),idx2(iyrs2),miss


!c   The first nrun-1 index values will be missing.
!c
      do 10 j = 1, nr-1
         idx(j) = miss
         idx2(j) = miss
 10   continue

!c
!c     Sum nrun precip. values; 
!c     store them in the appropriate index location.
!c
!c     If any value is missing; set the sum to missing.
!c
      do 30 j = nr, iyrs
         idx(j) = 0.0
          i=j-nr+1
          idx(j)=SUM(p(i:j))
 30   continue

      do 35 j = nr, iyrs2
         idx2(j) = 0.0
          i=j-nr+1
          idx2(j)=SUM(p2(i:j))
 35   continue
      return
      end


!!################################################# VARIANTE II ###################
!!Tomando toda la serie para construir la distribucion
! 
! 	ki=1
! 	do i=j,n,season   !delete the missing values.
!      if(dataSeries(i) .ne. rmiss)then
!       seasonSeries(ki)=dataSeries(i)
!      endif
!      ki=ki+1
!     enddo
! 	nSeason=ki-1
! 
!     k=nSeason+1
!     do i=j,nf,season   !delete the missing values.
!       if(dataSeries2(i).ne.rmiss)then
!       seasonSeries(k)=dataSeries2(i)
!     endif
!      k=k+1
!     enddo
! 	nSeason2=k-1
! 	
! 	
!     call ssort(seasonSeries, nSeason2)					!sort upward the data serie
! !     call pwm(seasonSeries, nSeason,beta, -0.35, 0.)
!     call pwm(seasonSeries, nSeason2,beta, 0., 0.)
!     call logLogisticFit(beta, logpara)
! 
!     do ik=1,3
!      logLogisticParams(j,ik)=logpara(ik)
!     enddo

!!############################################ FIN VARIANTE II ####################
! Calculates the first three probability weighted moments of a sample,
! using either the direct method from the definition of pwm (when A=B=0)
! or a plotting position formula (when A<=B<=0).

subroutine pwm(series,n,beta,A,B,isB)
	IMPLICIT NoNE
!   parameter (imaxyrs=250*12)
 integer i,inn,nn,n
 real acum(3), F,series(n),beta(3),A,B,rmiss,isB

 acum(1) = 0.
 acum(2) = 0.
 acum(3) = 0.
 if (A.eq.0..and.B.eq.0.)then
  do i=1,n
!    write(*,*) i,n,series(i)
    acum(1) = acum(1) +  series(i) 
    if(isB.eq.0.)then
     acum(2) = acum(2) + series(i) * (n-i)
     acum(3) = acum(3) + (series(i) * (n-i) * (n-i-1))
    endif
    if(isB.eq.1.)then
     acum(2) = acum(2) + series(i) * (i-1)
     acum(3) = acum(3) + (series(i) * (i-1) * (i-2))
    endif
  enddo
 else
  do i=1,n
    acum(1) = acum(1) +  series(i)
    F = (i+A) / (n+B)
    if(isB.eq.0.)then
     acum(2) = acum(2) + series(i) * (1-F)
     acum(3) = acum(3) + (series(i) * (1-F) * (1-F))
    endif
    if(isB.eq.1.)then
     acum(2) = acum(2) + series(i) * F
     acum(3) = acum(3) + (series(i) * F * F)
    endif
  enddo


!    F = (i-1+A) / (n+B)
!
!    acum(1) = acum(1)+series(i-1)
!    acum(2) = acum(2)+series(i-1)*(1-F)
!    acum(3) = acum(3)+series(i-1)*(1-F)*(1-F)
!  enddo
 endif
 beta(1) = acum(1) / n
 beta(2) = acum(2) / n / (n-1)    ! print*,params(2),value-params(1),v1,v3
 beta(3) = acum(3) / n / ((n-1)*(n-2))
 write(*,*)beta(1), beta(2), beta(3)
return
end

! logLogisticFit()
! Estimates the parameters of a Gamma distribution functions
subroutine logLogisticFit(beta,logLogisticParams)

 real g1, g2,beta(3),logLogisticParams(3)

! estimate gamma parameter
 logLogisticParams(3) = (2*beta(2)-beta(1)) / (6*beta(2)-beta(1)-6*beta(3))

 g1 = exp(gammLn1(1+1/logLogisticParams(3)))
 g2 = exp(gammLn1(1-1/logLogisticParams(3)))
! estimate alpha parameter
 logLogisticParams(2) = (beta(1)-2*beta(2))*logLogisticParams(3) / (g1*g2)
! estimate beta parameter
 logLogisticParams(1) = beta(1) - logLogisticParams(2)*g1*g2

return
end

! logLogisticCDF()
! Gives the cumulative distribution function of 'value', following a LogLogistic distribution

function logLogisticCDF(value, params)
 real value,params(3),v1,v2,v3,logLogisticCDF
 v1=params(2)/(value-params(1))
 v2=v1**params(3)
 v3=(1/(1+v2))
  ! write(*,*)'v1,v2,v3',v1,v2,v3
 logLogisticCDF=-standardGaussianInvCDF(v3)
!  write(*,*)'logLogisticCDF, logpara',logLogisticCDF, params
return
end

! gaussianInvCDF()
! Gives the inverse cumulative distribution function following a standard Gaussian
! distribution. I.e., the function gives the standardized value (Z) corresponding to
! a given probability (0<prob<1) following a standard Gaussian distribution (mean=0
! and sd=1).

function standardGaussianInvCDF(prob)
 real prob,C(3),d(4)
 real W, WW, WWW, results,standardGaussianInvCDF
 data C /2.515517,0.802853,0.010328/
 data d /0,1.432788,0.189269,0.001308/

 if (prob.le.0.5)then
   W = sqrt(-2*log(prob))
 else 
  W =sqrt(-2*log(1-prob))
 endif

 WW = W*W
 WWW = WW*W
 results = W - (C(1) + C(2)*W + C(3)*WW) / (1 + d(2)*W + d(3)*WW + d(4)*WWW)
 standardGaussianInvCDF=results
 if (prob.gt.0.5) standardGaussianInvCDF = -results
return 
end

!!c
!!c     For those who don't have a ln(gamma) function.
!!c
     function gammLn1(xx)
     dimension cof(6)
     data cof /76.18009173, -86.50532033, 24.01409822, -1.231739516,0.120858003e-2, -0.536382e-5/          
     x = xx - 1.0
     tmp = x + 5.5
     tmp = tmp - (x+0.5) * log (tmp)
     ser = 1.0
     do 10 j = 1, 6
        x = x + 1.0
        ser = ser + cof(j) / x
10   continue
     gammln = -tmp + log (2.50662827465 * ser/x)
     return
     end


	SUBROUTINE SSORT (X,N)
      IMPLICIT NONE
!
!    Example of a Sele!tion Sort   Using a Fortran 90 Intrinsi! Fun!tion
!
!***BEGIN PROLOGUE  SSORT
!***PURPOSE  Sort an array and make the same inter!hanges in
!            an auxiliary array.  The array is sorted in
!            de!reasing order.
!***TYPE      SINGLE PRE!ISION
!***KEYWORDS  SORT, SORTING
!
!   Des!ription of Parameters
!      X - array of values to be sorted   (usually abs!issas)
!      IY - array to be !arried with X (all swaps of X elements are
!          mat!hed in IY .  After the sort IY(J) !ontains the original
!          postition of the value X(J) in the unsorted X array.
!      N - number of values in array X to be sorted
!      KFLAG - Not used in this implementation
!
!***REVISION HISTORY  (YYMMDD)
!   950310  DATE WRITTEN
!   John Mahaffy
!***END PROLOGUE  SSORT
!     .. Scalar Arguments ..
      INTEGER KFLAG, N
!     .. Array Arguments ..  -----NOTE the 2 new ways of declaring array size
      REAL X(1:N)
      INTEGER IY(N)
!     .. Local Scalars ..
      REAL TEMP
      INTEGER I, ISWAP(1), ITEMP, ISWAP1
!     .. External Subroutines ..
!     None
!     .. Intrinsic Functions ..
      INTRINSIC MINLOC
!
!
!    MAXLOC is a FORTRAN 90 function that returns the index value for the
!    maximum element in the array
!***FIRST EXECUTABLE STATEMENT  SSORT
!
      DO 200 I=1,N-1

         ISWAP=MINLOC(X(I:N))

         ISWAP1=ISWAP(1)+I-1
         IF(ISWAP1.NE.I) THEN
           TEMP=X(I)
            X(I)=X(ISWAP1)
            X(ISWAP1)=TEMP
            ITEMP=IY(I)
            IY(I)=IY(ISWAP1)
            IY(ISWAP1)=ITEMP
         ENDIF
  200 CONTINUE
      RETURN
      END

! Sorts a given data series from the lowest to the highest value
subroutine upward(series,n)

 parameter (imaxyrs=250)
 integer i, j,n
 real series(imaxyrs),rmiss
 real*8 temp
 
 do i=0,n-1
  do j=0,n-1-i

    if (series(j+1).lt.series(j))then
     temp = series(j)
     series(j) = series(j+1)
     series(j+1) = temp
    endif

  enddo

 enddo

return
end


! Round value with idec decimals
function around(value, idec)
real value,ndec
integer idec
ndec=10.0**idec
around=FLOAT (INT(value * ndec + 0.5)) / ndec
return
end

