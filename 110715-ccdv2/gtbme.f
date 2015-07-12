c general 2body matrix elements
      subroutine gtbme(bast,leng,iwid,ns,nq,iobas,ham2)			!try !sijie

      implicit none

      integer :: leng,iwid
      integer :: i,j,k,l,m,a,b,bra,ket,ns,nq
      integer :: bast(leng,iwid),ham2(iwid,iwid)
      integer, dimension(ns,nq) :: iobas
      real    :: mel,melt


c i counters bras
c j counters ket
      do bra=1,iwid					
         do ket=1,iwid
            mel=0.d0
	      do i=1,leng
	        do j=1,leng
                do k=1,leng
                  do l=1,leng
                  if(bast(i,bra).ne.bast(j,bra).and.bast(k,ket)
     &.ne.bast(l,ket).and.iobas(bast(i,bra),1).eq.iobas(bast(j,bra),1)
     &.and.iobas(bast(k,ket),1).eq.iobas(bast(l,ket),1)
     &.and.iobas(bast(i,bra),2).eq.-1.and.iobas(bast(k,ket),2).eq.-1)
     &then
				melt=1.
	            else 
				melt=0.
				endif
                  a=1
                  b=1
c                  do m=1,2
c                    if(a.eq.i.or.a.eq.j)a=a+1
c                    if(b.eq.k.or.b.eq.l)b=b+1
c                    if(a.eq.i.or.a.eq.j)a=a+1
c                    if(b.eq.k.or.b.eq.l)b=b+1
c                    if(bast(a,bra).ne.bast(b,ket))melt=0.
c                    a=a+1
c                    b=b+1
c                  enddo
                  mel=mel+melt
                  enddo
                enddo
	        enddo
	      enddo
         ham2(bra,ket)=mel
         enddo
      enddo

      return
      end
