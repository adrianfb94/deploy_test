!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
!       Esta es una versión de calculo del SPI a partir de la versión 
!      desarrollada por la Universidad de Washington USA.
! 
! 
!        Se ha desarrollado un algoritmo para que sea ejecutada a partir de
!        una llamada combinada de scripts CSH-CDO-GrADS. Un Script  
!        GrADS crea los ficheros de entrada y el fichero de parametros 
!        que este programa utiliza para los calculos.
!        El file de entrada contiene el periodo base para el calculo de la distribución
!        y el periodo de datos que serán evaluados
! 
!        En esta versión se calcula el spi para el periodo de meses especificado en los parametros
!        Esta version es complementaria a la desarrollada por Roilan Hernandez
! 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
program spiprg
parameter (inlen=7)
character(len=:),allocatable::ifile,ifile2,ofile,result_file
character(len=:),allocatable::arg
real,allocatable,dimension(:,:,:)::tarray,spiout,pet,wtbalance,droughts
real,allocatable,dimension(:,:)::totumbral,conseumbral,longr,rndrought,spei,prcumbral
real beta(12),gamm(12),pzero(12),rumbral,sumbral
real miss,umbral,lat,dlat, mymaxvalue
real dlim(7)
real logLogisticParams,standardGaussianInvCDF
integer line(20),ix,iy,nyears,yr,inx,iny,i,k,ncd,tlen,tlen2,spilen,totlen,eta,etr,iarg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Lectura de los parametros pasados en linea de comando

!ichars=500! esto era para pasar a allocate de los nombre de files pero no funciono
allocate(character(len=500)::arg)
allocate(character (len=500) ::ifile,ifile2,ofile,result_file)! esto es para leer lasgos nombres de files
n=iargc() !n almacena en numero total de parametros en la cmdline

 call getarg(1,arg)
 ifile=arg !nombre fichero con datos de periodo de referencia
! write(*,*)ifile
 call getarg(2,arg)
 ifile2=arg !nombre fichero de datos de lluvia. Solo valido para el SPEI
! write(*,*)ifile2
 call getarg(3,arg)
 ofile=arg            !nombre de fichero salida del spi
! write(*,*)ofile
 call getarg(4,arg)
 result_file=arg      !nombre fichero con datos de conteo de sequias acorde al umbral
! write(*,*)result_file
do iarg=5,n-9
 call getarg(iarg,arg)
! write(*,*)iarg,arg
 read(arg,*) line(iarg-4) !parametros tlen,nmon,inx,iny,miss,tlen2,lenspi,ndrou
enddo
 call getarg(15,arg)
 read(arg,*) lat     !lat es la latitud inicial del fichero de datos de entrada
! write(*,*)lat
 call getarg(16,arg)
 read(arg,*) dlat    !dlat es el incremento de la latitud en el fichero
! write(*,*)dlat
 call getarg(17,arg)
 read(arg,*) line(11)
 call getarg(18,arg)
 read(arg,*) line(12)
 call getarg(19,arg)
 read(arg,*) line(13)
 call getarg(20,arg)
 read(arg,*) umbral     !este es el unico parametro que se lee fuera del ciclo por ser real
 call getarg(21,arg)
 read(arg,*) dlim(1)     
  call getarg(22,arg)
 read(arg,*) dlim(2)     
 call getarg(23,arg)
 read(arg,*) dlim(3)    
 dlim(4)=abs(dlim(3))
!*-1 
                       !si se quiere annadir mas parametros, se recomienda sea antes de este

!write(*,*) n,umbral,line(13)
!call sleep(5)

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


allocate(tarray(inx,iny,totlen)) ! allocate the array that should be read from input file

open(31,file=ifile,FORM='UNFORMATTED',access='direct',RECL=IRL)

IREC=1
do i=1,totlen
  read(31,REC=IREC) tarray(:,:,i)!read all the records (tarray(:,:,totlen)) as calculated by GrADS script
  IREC=IREC+1
enddo

ifirst=1
last=tlen2 !here tlen2 is the lenght of the assessment record

allocate(spiout(inx,iny,last),totumbral(inx,iny)) ! allocate spi and spiout variables
allocate(longr(inx,iny),rndrought(inx,iny),prcumbral(inx,iny))
allocate(droughts(inx,iny,inlen))
mymaxvalue=0.0

! Begging of the spi computation cycle
do iy=1,iny
 do ix=1,inx
   if(tlen.gt.tlen2)then
     nr=tlen
    else
     nr=tlen2
    endif
    if(tarray(ix,iy,itr).ne.miss.and.tarray(ix,iy,ita).ne.miss)then
     call spigam (spilen, nint(tarray(ix,iy,itr:etr))/100.0, nint(tarray(ix,iy,ita:eta))/100.0, &
                   beta, gamm, pzero, spiout(ix,iy,ifirst:last),miss,tlen,tlen2)
!   counting some total spi values below umbral
    totumbral(ix,iy)=0 !total of spi below umbral  
    longr(ix,iy)=0 ! maximum consecutive spi values below umbral
    rumbral=0 ! temporal counting variable for maximum consecutive spi values
    rndrought(ix,iy)=0
    sumbral=0
      do i=ifirst,last
         if(spiout(ix,iy,i)/=1.0e37 .and. mymaxvalue<=spiout(ix,iy,i))then
            mymaxvalue=spiout(ix,iy,i)
         end if
         ! if(spiout(ix,iy,i)==1.0e37)then
         !    write(*,*)ix,iy,i,spiout(ix,iy,i)
         !    ! spiout(ix,iy,i)=miss
         ! else
         !    if(mymaxvalue<=spiout(ix,iy,i))then
         !       mymaxvalue=spiout(ix,iy,i)
         !    end if
         ! end if
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
      prcumbral(ix,iy)=(totumbral(ix,iy)/(tlen2-spilen))*100
      do id=1,7
       droughts(ix,iy,id)=(droughts(ix,iy,id)/(tlen2-spilen))*100
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

deallocate(tarray)

do iy=1,iny
   do ix=1,inx
      do i=ifirst,last
         if(spiout(ix,iy,i)==1.0e37)then
            spiout(ix,iy,i)=mymaxvalue
         end if
      end do
   end do
end do


!c  Writing SPI values (spiout variable)

! write(*,*)minval(spiout), maxval(spiout)


open(49,file=ofile,access='direct',recl=IRL)
ir=1
do i=ifirst, last
  write(49,rec=ir) spiout(:,:,i) !se escriben los valores de spi para el area y periodo definido
  ir=ir+1
enddo
close(49)

! write(*,*) maxval(spiout),maxloc(spiout), mymaxvalue
! stop

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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c     These functions compute the Standardized Precipitation Index
!c     using an incomplete gamma distribution function to estimate
!c     probabilities.
!c
!c     Useful references are:
!c
!c     _Numerical Recipes in C_ by Flannery, Teukolsky and Vetterling
!c     Cambridge University Press, ISBN 0-521-35465-x 
!c
!c     _Handbook of Mathematical Functions_ by Abramowitz and Stegun
!c     Dover, Standard Book Number 486-61272-4
!c
!c     Notes for ForTran users:
!c     - There is a companion _Numerical Recipes in Fortran.  The
!c     following code was translated from the c code.
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c   Calculate indices assuming incomplete gamma distribution.
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spigam (nrun, pp, pp2, beta, gamm, pzero, index2,rmiss,maxyrs,maxyrs2)

      real pp(maxyrs), index(maxyrs),index2(maxyrs2),tmparr(maxyrs+1),beta(12), gamm(12), pzero(12),rmiss
      real pp2(maxyrs2)

      call cumsum (nrun,pp,pp2,index,index2,rmiss,maxyrs,maxyrs2)
!c
!c   For nrun<12, the monthly distributions will be substantially
!c   different.  So we need to compute gamma parameters for
!c   each month starting with the (nrun-1)th.
!c
      do 50 i = 0,11
         n = 0
         do 40 j = nrun+i, maxyrs, 12
            if(index(j) .ne. rmiss) then
               n = n + 1
               tmparr(n) = index(j)
            endif
 40      continue

         im = mod (nrun+i-1, 12) + 1

!c       
!c     Here's where we do the fitting.
!c
         call gamfit (tmparr,n,alpha,beta(im), gamm(im), pzero(im))

50   continue

!c
!c     Replace precip. sums stored in index with SPI's
!c
      do 60 j = nrun, maxyrs2
         im = mod (j-1,12) + 1

         
         if(index2(j) .ne. rmiss) then
!c
!c     Get the probability
!c
         
         index2(j) = gamcdf(beta(im), gamm(im), pzero(im), index2(j))
!c
!c     Convert prob. to z value. 
!c
         index2(j) = anvnrm(index2(j))
         ! if(index2(j)==1.0e37)then
         !    index2(j) = rmiss
         ! end if 
   
  
         endif
 60   continue
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


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c   input prob; return z.
!c
!c   See Abromowitz and Stegun _Handbook of Mathematical Functions_, p. 933
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function anvnrm (prob)
      data c0, c1, c2 /2.515517, 0.802853, 0.010328/
      data d1, d2, d3 /1.432788, 0.189269, 0.001308/
  
      if (prob .gt. 0.5) then
         sign = 1.0
         prob = 1.0 - prob
      else
         sign = -1.0
      endif
  
      if (prob .lt. 0.0) then
         anvnrm = 0.0

         return
      endif

      if (prob .eq. 0.0) then
         anvnrm = 1.0e37 * sign

         return
      endif
  
      t = sqrt(alog (1.0 / (prob * prob)))
      anvnrm = (sign * (t - ((((c2 * t) + c1) * t) + c0)/((((((d3 * t) + d2) * t) + d1) * t) + 1.0)))

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c  Estimate incomplete gamma parameters.
!c
!c  Input:
!c     datarr - data array
!c     n - size of datarr
!c
!c Output:
!c     alpha, beta, gamma - gamma paarameters
!c     pzero - probability of zero.
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gamfit (datarr, n, alpha, beta, gamm, pzero)
      real datarr(*)
      if (n .le. 0) then
         write(0, *) 'Error in gamfit - empty data array'
         stop
      endif

      sum = 0.0
      sumlog = 0.0
      pzero = 0.0
      nact = 0

!c     compute sums
      do 10 i = 1, n
         if (datarr(i) .gt. 0.0) then
            sum = sum + datarr(i)
            sumlog = sumlog + alog (datarr(i))
            nact = nact + 1
         else
            pzero = pzero + 1
         endif
 10   continue
      pzero = pzero / n
      if(nact .ne. 0.0) av = sum / nact
  
!c     Bogus data array but do something reasonable
      if(nact .eq. 1) then
         alpha = 0.0
         gamm = 1.0
         beta = av
         return
      endif

!c     They were all zeroes. 
      if(pzero .eq. 1.0) then
         alpha = 0.0
         gamm = 1.0
         beta = av
         return
      endif

!c     Use MLE
      alpha = alog (av) - sumlog / nact
      gamm = (1.0 + sqrt (1.0 + 4.0 * alpha / 3.0)) / (4.0 * alpha)
      beta = av / gamm
  
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c  Compute probability of a<=x using incomplete gamma parameters.
!c
!c  Input:
!c      beta, gamma - gamma parameters
!c      pzero - probability of zero.
!c      x - value.
!c
!c  Return:
!c      Probability  a<=x.
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function  gamcdf (beta, gamm, pzero, x)
      if(x .le. 0.0) then
         gamcdf = pzero
      else
!        write(*,*) gammap(gamm, x / beta)
         gamcdf = pzero + (1.0 - pzero) * gammap (gamm, x / beta)
      endif
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c  Compute inverse gamma function i.e. return x given p where CDF(x) = p.
!c
!c  Input:
!c      beta, gamma - gamma parameters
!c      pzero - probability of zero.
!c      prob - probability.
!c
!c  Return:
!c      x as above.
!c
!c  Method:
!c      We use a simple binary search to first bracket out initial
!c      guess and then to refine our guesses until two guesses are within
!c      tolerance (eps).  Is there a better way to do this?
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function gaminv (beta, gamm, pzero, prob)
      data  eps /1.0e-7/
  
!c     Check if prob < prob of zero
      if (prob .le. pzero) then
         gaminv = 0.0
         return
      endif
  
!c     Otherwise adjust prob
      prob = (prob - pzero) / (1.0 - pzero)
  
!c     Make initial guess. Keep doubling estimate until prob is
!c     bracketed.
      thigh = 2.0*eps 
 10   continue
      phigh = gamcdf (beta, gamm, pzero, thigh)
      if(phigh .ge. prob) goto 20
      thigh = thigh*2.0
      goto 10
 20   continue
      tlow = thigh / 2.0
  
!c     Iterate to find root.
      niter = 0
 30   continue
      if((thigh - tlow) .le. eps) goto 40
      niter = niter + 1
      t = (tlow + thigh) / 2.0
      p = gamcdf (beta, gamm, pzero, t)
      
      if (p .lt. prob) then
         tlow = t
      else
         thigh = t
      endif
      goto 30
 40   continue
      gaminv = (tlow + thigh) / 2.0
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c  Functions for the incomplete gamma functions P and Q
!c
!c                  1     /x  -t a-1
!c   P (a, x) = -------- |   e  t    dt,  a > 0
!c              Gamma(x)/ 0
!c
!c   Q (a, x) = 1 - P (a, x)
!c
!c Reference: Press, Flannery, Teukolsky, and Vetterling, 
!c        _Numerical Recipes_, pp. 160-163
!c
!c Thanks to kenny@cs.uiuc.edu
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c Evaluate P(a,x) by its series representation.  
!c
      function gamser (a, x)
!c     Maximum number of iterations, and bound on error.
      parameter (maxitr=100, eps=3.0e-7)
      data iwarn /0/
      gln = gammln (a)
      if (x .eq. 0.0) then
         gamser = 0.0
         return
      endif
      ap = a
      sum = 1.0 / a
      del = sum
  
      do 10 n = 1, maxitr
         ap = ap + 1.0
         del = del * (x / ap)
         sum = sum + del
         if (abs (del) .lt. eps * abs (sum)) goto 20
 10   continue
      iwarn = iwarn + 1
      if (iwarn .lt. 20) then
!         write (0, *) 'gamser(',a,x,'): not converging.'
!         write (0, *) 'Approximate value of ',sum,'  + /-',del,' used.'
      endif
 20   continue
      gamser =  sum * exp (-x + a * alog (x) - gln)
      return
      end

!c
!c     Evaluate P(a,x) in its continued fraction representation.
!c
      function gammcf (a, x)
      parameter (maxitr=100, eps=3.0e-7)
      data nwarn / 0 /, g / 0.0 /
  
      gln = gammln (a)
      gold = 0.0
      a0 = 1.0
      a1 = x
      b0 = 0.0
      b1 = 1.0
      fac = 1.0
      do 10 n = 1, maxitr
         an = n
         ana = an - a
         a0 = (a1 + a0 * ana) * fac
         b0 = (b1 + b0 * ana) * fac
         anf = an * fac
         a1 = x * a0 + anf * a1
         b1 = x * b0 + anf * b1
         if (a1 .ne. 0.0) then
            fac = 1.0 / a1
            g = b1 * fac
            if (abs((g - gold) / g) .lt. eps) goto 20
            gold = g
         endif
 10   continue
      nwarn = nwarn + 1
      if (nwarn .lt. 20) then
      endif
 20   continue

      gammcf =  g * exp (-x + a * alog (x) - gln)

      return
      end
!c
!c     Evaluate the incomplete gamma function P(a,x), choosing the most 
!c     appropriate representation.
!c
      function gammap (a, x)
      if (x .lt. a + 1.0) then
         gammap = gamser (a, x)
      else
         gammap = 1.0 - gammcf (a, x)
      endif
      return
      end
!c
!c     Evaluate the incomplete gamma function Q(a,x), choosing the most 
!c   appropriate representation.
!c
      function gammaq (a, x)
      if (x .lt. a + 1.0) then
         gammaq = 1.0 - gamser (a, x)
      else
         gammaq = gammcf (a, x) 
      endif
      return
      end
!c
!c     For those who don't have a ln(gamma) function.
!c
      function gammln(xx)
      dimension cof(6)
      data cof /76.18009173, -86.50532033, 24.01409822, -1.231739516,0.120858003e-2, -0.536382e-5/          
      x = xx - 1.0
      tmp = x + 5.5
      tmp = tmp - (x+0.5) * alog (tmp)
      ser = 1.0
      do 10 j = 1, 5
         x = x + 1.0
         ser = ser + cof(j) / x
 10   continue
      gammln = -tmp + alog (2.50662827465 * ser)
      return
      end
