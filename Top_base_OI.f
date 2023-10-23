      real*8 function sigox(il,e1)
c     interpolates the O I cross sections from TopBase for the different
c      levels levi
c     Note LS coupling so that FS levels 1-3 are all in level 1
c     cross setion in MBarns
c     CF 230202      
      implicit none
      integer nltb,nlcr,ncr,inittb,levi,il,id
      common/inttopbase/inittb
      parameter(nltb=500,nlcr=1000,ncr=1000)
      integer i,j,k,imax,n,npi(nltb),nmin,lev(nltb),srr(nltb),np(nltb),
     &     ll,g,ne
      real*8 e00(nltb),en(500,nlcr),cross(500,nlcr),e1,cr,enry,
     &     slope(500)      
      character dum
      save en,cross,npi,slope
c correct for first level in LS is a triplet. rest shifted up to at least 13
      if(il <= 3) then
         levi=1
      else
         levi=il-2
      endif
      
      if(inittb==1) then
         inittb=0
         open(11,file='./ATDAT/O_I_TopBase_lev_lines.dat',status='old')
         do i=1,1000
            read(11,*,err=11,end=11)srr(i),lev(i),ll,ll,ll,id,g,e00(i),
     &           np(i)
         enddo
 11      imax=i-1
         close(11)
         open(11,file='./ATDAT/O_I_TopBase.dat',status='old')
         do k=1,15
            do i=1,imax
               if(k==lev(i)) then
                  nmin=srr(i)
                  npi(k)=np(i)
               endif
            enddo
            if(k==1) then
               do n=1,4
                  read(11,*)dum
               enddo
            else
               read(11,*)dum
            endif
            do j=1,npi(k)
               read(11,*)en(k,j),cross(k,j)
               enry=en(k,j)
               en(k,j)=13.6057*en(k,j)
            enddo
            slope(k)=log10( cross(k,npi(k)) / cross(k,npi(k)-10) ) /
     &           log10( en(k,npi(k)) / en(k,npi(k)-10) )
         enddo
      endif
      cr=0.
      if(e1 > en(levi,npi(levi)) ) then
         cr=cross(levi,npi(levi))*(e1/en(levi,npi(levi)))**slope(levi)
      else
         do j=1,npi(levi)
            if (e1 > en(levi,j) .and. e1 <= en(levi,j+1)) then
               cr = (cross(levi,j) + cross(levi,j+1))/2.
               goto 33
            endif
         enddo
 33      continue
      endif
      sigox=cr
      return
      end
      
