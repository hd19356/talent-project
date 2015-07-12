      call lccd
      stop
      end

      subroutine lccd

      implicit none
c______
      integer :: i,idgn,ig,j,k,l,m,n2p,n2h,nl,nq,ns
      integer :: iwid,iwid2,iepsp,iepsm,ind
c______
      real    :: dg,g,g2,deccd
c______
      parameter(ns=8,nq=3)
c______
      real,    dimension(4)     :: eps,b,x
c______
      integer, dimension(ns,nq) :: iobas
      integer, dimension(2,28)  :: itbas,ic2p,ic2h
      integer, dimension(2,2)   :: ic2pp,ic2hp
      integer, dimension(2,4)   :: iph2b
      integer, dimension(4,4)   :: ham2
c______
      real,    dimension(4,4)   :: vmat,aa
      real,    dimension(4,5)   :: alccd
c______
      dg   =0.01
      g    =-1.

 1000 format(50i4)
 1001 format(50f8.3)
 1002 format('interation for g=',f8.3)

      open(unit=1,file='lccd.inp',status='old')
      open(unit=2,file='lccd.out',status='unknown')

      read(1,*)nl,idgn
      close(1)

c # of states = # of energy levels X degeneracy
      write(*,*)
      write(*,*)' # of levels, degeneracy, total # of states:'
      write(*,*)
      write(*,*)nl,idgn,ns

c build 1body basis
      call obas(nl,ns,nq,iobas)

c build 2body basis
      call tbas(ns,nq,iobas,itbas)

c store 2p and 2h configurations
      ic2p(:,:)=0
      ic2h(:,:)=0
      k=0
      l=0

      do i=1,1
       do j=1,28
        if(itbas(i,j).le.ns/2.and.itbas(i+1,j).le.ns/2) then
         ic2h(i,j)=itbas(i,j)
         ic2h(i+1,j)=itbas(i+1,j)
         if(abs(ic2h(i+1,j)-ic2h(i,j)).eq.1) then
          if(ic2h(i,j)/2*2.ne.ic2h(i,j)) k=k+1
         endif 
        endif
        if(itbas(i,j).gt.ns/2.and.itbas(i+1,j).gt.ns/2) then
         ic2p(i,j)=itbas(i,j)
         ic2p(i+1,j)=itbas(i+1,j)
         if(abs(ic2p(i+1,j)-ic2p(i,j)).eq.1) then
          if(ic2p(i,j)/2*2.ne.ic2p(i,j)) l=l+1
         endif
        endif
       enddo
      enddo
      n2p=k
      n2h=l
 
      write(*,*)
      write(*,*)'2h configurations:'
      write(*,*)
      do i=1,2
       write(*,1000)(ic2h(i,j),j=1,28)
      enddo
      write(*,*)
      write(*,*)('2p configurations:')
      write(*,*)
      do i=1,2
       write(*,1000)(ic2p(i,j),j=1,28)
      enddo
      write(*,*)

c retain only paired configurations
      k=1
      l=1
      do i=1,1
       do j=1,28
        if(ic2p(i,j).ne.0) then
         if(abs(ic2p(i+1,j)-ic2p(i,j)).eq.1) then
          if(ic2p(i,j)/2*2.ne.ic2p(i,j)) then
           ic2pp(i,k)=ic2p(i,j)
           ic2pp(i+1,k)=ic2p(i+1,j)
           k=k+1
          endif
         endif
        endif
        if(ic2h(i,j).ne.0) then
         if(abs(ic2h(i+1,j)-ic2h(i,j)).eq.1) then
          if(ic2h(i,j)/2*2.ne.ic2h(i,j)) then
           ic2hp(i,l)=ic2h(i,j)
           ic2hp(i+1,l)=ic2h(i+1,j)
           l=l+1
          endif
         endif
        endif
       enddo
      enddo

      write(*,*)
      write(*,*)'2-hole paired configurations:'
      write(*,*)
      do i=1,2
       write(*,1000)(ic2hp(i,j),j=1,n2p)
      enddo
      write(*,*)
      write(*,*)'2-particle paired configurations:'
      write(*,*)
      do i=1,2
       write(*,1000)(ic2pp(i,j),j=1,n2h)
      enddo
      write(*,*)

c build interaction matrix
      iph2b(:,:)=0

      do i=1,2
         do j=1,n2p
            iph2b(i,j)=ic2hp(i,j)
         enddo
         do j=1,n2h
            iph2b(i,j+n2p)=ic2pp(i,j)
         enddo
      enddo

      write(*,*)
      write(*,*)('paired basis:')
      write(*,*)
      do i=1,2
         write(*,1000)(iph2b(i,j),j=1,n2p+n2h)
      enddo
      write(*,*)

      iwid =n2p+n2h
      iwid2=iwid/2     

      call gtbme(iph2b,2,iwid,ns,nq,iobas,ham2)

      write(*,*)
      write(*,*)'interaction matrix:'
      write(*,*)
      do i=1,iwid
         write(*,1000)(ham2(i,j),j=1,iwid)
      enddo
      write(*,*)

c set up lccd matrix only 2p and 2h contributions
      alccd(:,:)=0

c set up epsilon vector for lccd matrix
      i=1
      k=1
      l=1
      do m=1,iwid2
         do j=1,iwid2
            iepsp  =iobas(ic2pp(i,l),1)+iobas(ic2pp(i+1,l),1)
            iepsm  =iobas(ic2hp(i,j),1)+iobas(ic2hp(i+1,j),1)
            eps(k) =1.*(iepsp-iepsm)
            k=k+1
         enddo
         l=l+1
      enddo
c interaction strength loop
      do ig=1,200
         g =g+dg
         g2=g/2.
         deccd=0.
         write(*,*)
         write(*,1002)g
         write(*,*)

         do i=1,iwid
            do j=1,iwid
               vmat(i,j)=-1.*g2*ham2(i,j)
            enddo
         enddo

         write(*,*)
         write(*,*)'full V matrix:'
         write(*,*)
         do i=1,iwid
            write(*,1001)(vmat(i,j),j=1,iwid)
         enddo
         write(*,*)

c add interaction matrix elements to epsilon vector
         i=1
         j=1
         do k=1,4
            if(k.lt.iwid2+1) then 
               eps(k)=eps(k)+vmat(iwid2+1,iwid2+1)
               eps(k)=eps(k)+vmat(i,i)
               i=i+1
            endif
            if(k.ge.iwid2+1) then 
               eps(k)=eps(k)+vmat(iwid,iwid)
               eps(k)=eps(k)+vmat(j,j)
               j=j+1
            endif
         enddo   

c build lccd matrix
         k=3
         do i=1,iwid,2
            alccd(i,1)=vmat(k,l)
            alccd(i+1,1)=vmat(k,l+1)
            k=k+1
         enddo

         do i=1,iwid
            alccd(i,i+1)=eps(i)
         enddo

         j=iwid2+1
         k=iwid+1
         do i=1,iwid
            if(i.le.iwid2) then 
               alccd(i,i+3)=vmat(iwid2+1,iwid)
               alccd(i,j)=vmat(j-1,i)
               j=j-1
            endif
            if(i.gt.iwid2) then 
               alccd(i,i-1)=vmat(iwid,iwid-1)
               alccd(i,k)=vmat(k-3,i-2)
               k=k-1
            endif
         enddo         

         write(*,*)
         write(*,*)'lccd matrix:'
         write(*,*)

         do i=1,iwid
            write(*,1001)(alccd(i,j),j=1,iwid+1)	
         enddo

c set up non-homogenous system
         x(:)=0
         do i=1,iwid
            b(i)=-1.*alccd(i,1)
            do j=1,iwid
               aa(i,j)=alccd(i,j+1)
            enddo
         enddo

         write(*,*)
         write(*,*)'square matrix:'
         write(*,*)
         do i=1,iwid
            write(*,1001)(aa(i,j),j=1,iwid)
         enddo

         write(*,*)
         write(*,*)'non-homogenous term:'
         write(*,*)
         do i=1,iwid
            write(*,1001)b(i)
         enddo
         write(*,*)

c solve non-homogenous system
         call gauss(iwid,iwid,aa,b,x,ind)

         write(*,*)
         write(*,*)'solution of the non-homogenous system:'
         write(*,*)
         do i=1,iwid
            write(*,1001)x(i)
         enddo
         write(*,*)

c calculate correlation energy
         do i=1,iwid
            if(i.le.2) deccd=deccd+0.25*vmat(i,3)*x(i)
            if(i.gt.2) deccd=deccd+0.25*vmat(i,4)*x(i)
         enddo
         write(2,1001)g,deccd,x
c         write(*,1001)g,deccd,x
      enddo
      write(*,*)

      close(2)

      return
      end
