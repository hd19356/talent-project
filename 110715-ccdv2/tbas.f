      subroutine tbas(ns,nq,iobas,itbas)

      implicit none
      integer :: i,isum,j,k,l,m,np,ns,nsp,nq,numc
      integer, dimension(2,28)   :: itbas
      integer, dimension(ns)    :: i2bs
      integer, dimension(ns,nq) :: iobas
      isum=0
      k   =0
      l   =0
      m   =0
      np  =2
      nsp=ns-np
      numc=1

 1000 format('wrong number of particles:',i4)
 1001 format('wrong number of configurations:',2i4)
 2000 format(50i4)

c      write(*,*)'test'
c      i2bs(:)=0

      write(*,*)
      write(*,*)'all 2body configurations:'
      write(*,*)
      do i=1,ns
c         write(*,*)i
         do j=i+1,ns
            i2bs(:)=0
            i2bs(i)=1
            i2bs(j)=1
            do k=1,ns
               isum=isum+i2bs(k)
            enddo
            if(isum.ne.2) then
               write(*,1000)isum
               stop
            endif
            write(*,*)i2bs
            l=l+1
            isum=0
            if(iobas(i,1).eq.iobas(j,1)) then
               if(iobas(i,2).eq.iobas(j,2)) goto 10
            endif
            m=m+1
            itbas(1,m)=i
            itbas(2,m)=j
c            write(*,*)itbas(2,m)
   10       continue
         enddo
      enddo


      if(l.ne.28) then
         write(*,*)
         write(*,1001)l,28
         write(*,*)
         stop
      endif

      write(*,*)
      write(*,*)'2b-basis:'
      write(*,*)

      do i=1,2
         write(*,2000)(itbas(i,j),j=1,28)
      enddo
      write(*,*)

      return
      end
