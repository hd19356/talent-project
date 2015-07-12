      subroutine obas(nl,ns,nq,iobas)

      implicit none
      integer :: i,j,nl,nq,ns
      integer, dimension(ns,nq) :: iobas
      i=0

 1000 format('wrong number of states:',2i4)

c all 1b quantum numbers
      do j=1,nl
         i=i+1
         iobas(i,1)=j-1
         iobas(i,2)=1
         iobas(i,3)=i
         i=i+1
         iobas(i,1)=j-1
         iobas(i,2)=-1
         iobas(i,3)=i
      enddo

      if(i.ne.ns) then
         write(*,1000)i,ns
         stop
      endif
  
      write(*,*)
      write(*,*)'1body states:'
      write(*,*)
      do i=1,ns
         write(*,*)(iobas(i,j),j=1,3)
      enddo
      write(*,*)

      return
      end
