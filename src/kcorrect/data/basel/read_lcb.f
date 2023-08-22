      program read_lcb

c     -----------------------------------------------------------------------
c     Read wavelengths, model fluxes and model parameters in a file designed by
c     the metallicity, [Fe/H]
c     parameters   : Teff, logg and [Fe/H]     
c     wave         : 1221 wavelength points     in nm
c     Hnu          : 1221 flux moment points per model in erg/cm2/s/nm/sr
c     -----------------------------------------------------------------------

      implicit none

      integer 	      i,j,im
      integer 	      nmod,nmodel,iteff
      real*4          glog,feh,mhx
      real*8          wave(1221),hnu(1221,600)

      integer         nmet
      parameter       (nmet=19)

      real*4          met(nmet)
      data met        /-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,
     &                 -0.5,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.5,+0.1/

      character*3     mcod(nmet),suff,type
      data mcod       /'m50','m45','m40','m35','m30','m25','m20','m15',
     &                 'm10','m05','m03','m02','m01','p00','p01','p02',
     &                 'p03','p05','p10'/

      character*80    lfile  

c     ----------------------------------------------------------------  

c     assign metallicity file
c     -----------------------

      write(6,'(x,a,$)') 'Enter metallicity : '
      read(5,*)mhx    
      write(6,'(x,a,$)') 'Original/Corrected spectra ? (cor/ori) : '
      read(5,'(a3)')type

      do im = 1,nmet
         if(mhx.eq.met(im)) suff = mcod(im)
      enddo

      lfile = 'lcb'//suff//'.'//type

c     read models
c     -----------

      open(unit=1,file=lfile,status='old')
         
      READ(1,10) (WAVE(i),i=1,1221)
         
      DO 1000 j=1,600
         READ(1,11,end=999)NMOD,ITEFF,GLOG,FEH
         READ(1,12)(Hnu(i,j),i=1,1221)
 1000 CONTINUE

 999  close(unit=1)
      nmodel = j
         
         
 10   FORMAT(8f10.1)
 11   FORMAT(I6,I6,2f6.2)
 12   FORMAT(7e11.4)
      
      END
