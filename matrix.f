      subroutine cons(nd,n,c,a)
c
c move constant c ---> a
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
      do i=1,n
        do j=1,n
          a(i,j)=c
        enddo
      enddo
      return
      end

      subroutine zero(nd,n,a)
c
c 0 ---> a
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
      call cons(nd,n,0.d0,a)
      return
      end

      subroutine move(nd,n,a,b)
c
c move matrix a ---> b
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd)
      do i=1,n
        do j=1,n
          b(i,j)=a(i,j)
        enddo
      enddo
      return
      end

      subroutine movev(n,a,b)
c
c move vector a ---> b
c
      implicit double precision (a-h,o-z)
      dimension a(n),b(n)
      do i=1,n
        b(i)=a(i)
      enddo
      return
      end

      subroutine sum(nd,n,a,b,c)
c
c matrix summation c = a + b
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd),c(nd,nd)
      do i=1,n
        do j=1,n
          c(i,j)=a(i,j)+b(i,j)
        enddo
      enddo
      return
      end

      subroutine dif(nd,n,a,b,c)
c
c matrix difference c = a - b
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd),c(nd,nd)
      do i=1,n
        do j=1,n
          c(i,j)=a(i,j)-b(i,j)
        enddo
      enddo
      return
      end

      subroutine prod(nd,n,a,b,c)
c
c matrix product c = a * b
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd),c(nd,nd)
      do i=1,n
        do j=1,n
          c(i,j)=0.
          do k=1,n
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
          enddo
        enddo
      enddo
      return
      end

      subroutine prov(nd,n,a,b,c)
c
c multiply matrix by a vector c = a * b
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd),c(nd)
      do i=1,n
        c(i)=0.
        do j=1,n
          c(i)=c(i)+a(i,j)*b(j)
        enddo
      enddo
      return
      end

      subroutine mult(nd,n,c,a,b)
c
c multiply the constant c by a matrix a
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd)
      do i=1,n
        do j=1,n
          b(i,j)=c*a(i,j)
        enddo
      enddo
      return
      end

      subroutine trans(nd,n,a,b)
c
c transpose a matrix
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd)
      do i=1,n
        do j=1,n
          b(i,j)=a(j,i)
        enddo
      enddo
      return
      end

      subroutine primat(nf,nd,n,a)
c
c print a matrix in the file "nf"
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
 2000 format(' ')
 2010 format(20e12.4)
      write(nf,2000)
      do i=1,n
        write(nf,2010)(a(i,j),j=1,n)
      enddo
      return
      end

      subroutine delmat(nd,n,a)
c
c move delta function into a matrix
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
      do i=1,n
        do j=1,n
          a(i,j)=deltij(i,j)
        enddo
      enddo
      return
      end

      function vmax(nd,n,a)
c
c maximal absolute value in the matrix a
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
      vmax=0.
      do i=1,n
        do j=1,n
          if(abs(a(i,j)).gt.abs(vmax))vmax=a(i,j)
        enddo
      enddo
      return
      end

      function vmin(nd,n,a)
c
c minimal absolute value in the matrix a
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
      vmin=0.
      do i=1,n
        do j=1,n
          if(abs(a(i,j)).lt.abs(vmin))vmin=a(i,j)
        enddo
      enddo
      return
      end

      subroutine compr(nd,n,a,b)
c
c compress matrix a into b
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(1)
      do j=1,n
        do i=1,n
          k=j*(n-1)+i
          b(k)=a(i,j)
        enddo
      enddo
      return
      end

      subroutine decom(nd,n,b,a)
c
c decompress matrix b into a
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(1)
      do j=1,n
        do i=1,n
          k=j*(n-1)+i
          a(i,j)=b(k)
        enddo
      enddo
      return
      end

      subroutine double(nd,n,a,b)
c
c multiplexing the matrix a into b
c order of matrices is n*2
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd)
      do i=1,n
        in=i+n
        i1=(i-1)*2+1
        i2=i1+1
        do j=1,n
          jn=j+n
          j1=(j-1)*2+1
          j2=j1+1
          b(i1,j1)=a(i ,j )
          b(i1,j2)=a(i ,jn)
          b(i2,j1)=a(in,j )
          b(i2,j2)=a(in,jn)
        enddo
      enddo
      return
      end

      subroutine dedoub(nd,n,b,a)
c
c demultiplexing the matrix b into a
c order of matrices is n*2
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd)
      do i=1,n
        in=i+n
        i1=(i-1)*2+1
        i2=i1+1
        do j=1,n
          jn=j+n
          j1=(j-1)*2+1
          j2=j1+1
          a(i ,j )=b(i1,j1)
          a(i ,jn)=b(i1,j2)
          a(in,j )=b(i2,j1)
          a(in,jn)=b(i2,j2)
        enddo
      enddo
      return
      end

      subroutine tesdel(nd,n,a,ind)
c
c test a matrix to be equal with the delta function
c ind   =  0 : positive test
c     .ne. 0 : negative test
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
 9000 format(' *** a(',i3,',',i3,') =',d12.4)
      small=1.d-2
      ind=0
      do i=1,n
        do j=1,n
        if(dabs(a(i,j)-deltij(i,j)).gt.small) then
          ind=ind+1
c          write(* ,9000)i,j,a(i,j)
c          write(20,9000)i,j,a(i,j)
        endif
        enddo
      enddo
      return
      end

      subroutine tesdel2(nd,n,a,ind)
c
c test a matrix to be equal with the 2x2 block diagonal form
c ind = 0 : yes
c       n : no
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
      small=1.d-5
      ind=0
      do i=1,n
        do j=1,n
          if(iabs(i-j).gt.1) then
            if(dabs(a(i,j)).gt.small) ind=ind+1
          endif
        enddo
      enddo
      return
      end

      subroutine tesnul(nd,n,a,ind)
c
c test columns and rows of a matrix to be equal with 0
c ind = 0 : yes
c       n : no
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
      small=1.d-20
c
c test columns
      do i=1,n
        ind=0
        do j=1,n
          if(dabs(a(j,i)).gt.small) ind=ind+1
        enddo
        if(ind.eq.0) return
      enddo
c
c test rows
      do i=1,n
        ind=0
        do j=1,n
          if(dabs(a(i,j)).gt.small) ind=ind+1
        enddo
        if(ind.eq.0) return
      enddo
      return
      end

      function deltij(i,j)
      implicit double precision (a-h,o-z)
      deltij=0.d0
      if(i.eq.j) deltij=1.d0
      return
      end

      subroutine gauss(nd,n,a,b,x,ind)
c
c solve a system of linear equations using Gauss method
c
c nd  : dimension
c n   : order of the matrix
c a   : matrix of the coefficients
c b   : free terms vector
c x   : solution vector
c ind = 0 : compatible system
c       1 : incompatible system
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd),x(nd)
      ind=0
      im=n-1
      do i=1,im
        m1=i
        do m2=i,n
          if(abs(a(m1,i)).lt.abs(a(m2,i))) m1=m2
        end do
        if(a(m1,i).eq.0) then
          ind=1
          go to 90
        end if
        if(m1.eq.i) go to 25
        do j=1,n
          y=a(i,j)
          a(i,j)=a(m1,j)
          a(m1,j)=y
        end do
        y=b(i)
        b(i)=b(m1)
        b(m1)=y
   25   y=a(i,i)
        if(y.eq.0) then
          ind=1
          go to 90
        end if
        do j=i,n
          a(i,j)=a(i,j)/y
        end do
        b(i)=b(i)/y
        jm=i+1
        do j=jm,n
          if(a(j,i).eq.0) go to 40
          y=a(j,i)
          do k=i,n
            a(j,k)=a(j,k)-a(i,k)*y
          end do
          b(j)=b(j)-b(i)*y
   40     continue
        end do
      end do
      if(a(n,n).eq.0) then
        ind=1
        go to 90
      end if
      x(n)=b(n)/a(n,n)
      nm=n-1
      do i=1,nm
        x(n-i)=b(n-i)
        do k=1,i
          x(n-i)=x(n-i)-a(n-i,n-i+k)*x(n-i+k)
        end do
      end do
   90 continue
      return
      end

      subroutine invmat (nd,n,a,b,ind,det,idet)
c
c invert a matrix using Gauss method
c nd  : dimension
c n   : order of the matrix
c a   : matrix
c b   : inverse matrix
c ind : 0 : no inversion
c det : determinat of a
c idet: 0 : does not compute the determinant
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd),l(500)
      do i=1,n
        do j=1,n
          b(i,j)=0.
        end do
      end do
      do i=1,n
        b(i,i)=1.
        l(i)=i
      end do
      det=1.
      do i=1,n
        if(a(i,i).eq.0.) then
    8     k=i+1
    9     if(k.gt.n) go to 11
          if(a(i,k).ne.0) go to 12
          k=k+1
          go to 9
   11     ind=0
          return
   12     ia=l(i)
          l(i)=l(k)
          l(k)=ia
          do j=1,n
            aa=b(i,j)
            b(i,j)=b(k,j)
            b(k,j)=aa
            aa=a(i,j)
            a(i,j)=a(k,j)
            a(k,j)=aa
          end do
          go to 3
        endif
    3   x=a(i,i)
        if(idet.ne.0)det=det*x
c        if(dabs(det).gt.1.d+100)det=det/1.d+100
c        if(dabs(det).lt.1.d-100)det=det*1.d+100
        do k=1,n
          a(i,k)=a(i,k)/x
          b(i,k)=b(i,k)/x
        end do
        do k=1,n
          if(k.eq.i) go to 6
          x=a(k,i)
          if(x.eq.0) go to 6
          do j=1,n
            a(k,j)=a(k,j)-a(i,j)*x
            b(k,j)=b(k,j)-b(i,j)*x
          end do
    6     continue
        end do
      end do
      do i=1,n
        do j=1,n
          a(l(i),j)=b(i,j)
        end do
      end do
      ind=1
      return
      end

      subroutine homog(nd,n,k,a,b,x,ind)
c
c solve a homogeneous system of equations : a*x = 0
c k = used line to solve the system
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd),b(nd,nd),x(nd),c(500)
 9000 format(' *** ind =',i5)
      n1=n-1
      if(n1.le.0)return
c     if(k.eq.n) go to 10
      if(k.gt.n) go to 10
      do i=1,n
        c(i)=-a(i,k)
      enddo
      do i=1,n
        a(i,k)=0.d0
        a(k,i)=0.d0
      enddo
      a(k,k)=1.d0
      c(k)=1.d0
      call invmat(nd,n,a,b,ind,det,1)
      if(ind.eq.0) then
        write(*,9000)
        return
      endif
      call prov(nd,n,b,c,x)
      go to 20
   10 continue
      do i=1,n1
        c(i)=-a(i,n)
      enddo
      call invmat(nd,n1,a,b,ind,det,1)
      if(ind.eq.0) then
        write(*,9000)
        return
      endif
      call prov(nd,n1,b,c,x)
      x(n)=1.d0
c
c normalize solution
   20 continue
      s=0.d0
      do i=1,n
        s=s+x(i)**2
      enddo
      s=dsqrt(s)
      do i=1,n
        x(i)=x(i)/s
      enddo
      return
      end

      subroutine invers(nd,n,a,ainv,adjug,det,p)
c
c invert matrix a(nd,nd) of dimension n
c ainv  = inverse matrix
c adjug = det*ainv = adjoint matrix
c det   = determinat
c p(nd) = coefficients of the characteristic polinomial
c
      implicit double precision(a-h,o-z)
      dimension a(nd,nd),ainv(nd,nd),adjug(nd,nd),p(nd)
      small=1.d-30
      call move(nd,n,a,ainv)
      do k=1,n
        p(k)=0.d0
        do i=1,n
          p(k)=p(k)+ainv(i,i)
        enddo
        p(k)=p(k)/dfloat(k)
        if(k.eq.n) go to 5
        call move(nd,n,ainv,adjug)
        do i=1,n
          adjug(i,i)=ainv(i,i)-p(k)
        enddo
        call prod(nd,n,a,adjug,ainv)
      enddo
    5 continue
      call move(nd,n,adjug,ainv)
      if(dabs(p(n)).ge.small) then
        coef=1.d0/p(n)
        call scaled(nd,n,coef,ainv)
      endif
      det=p(n)
      if(mod(n,2).eq.1) return
      det=-det
      do i=1,n
        do j=1,n
          adjug(i,j)=-adjug(i,j)
        enddo
      enddo
      return
      end

      subroutine scaled(nd,n,c,a)
c
c a ---> c*a
c
      implicit double precision (a-h,o-z)
      dimension a(nd,nd)
      do i=1,n
        do j=1,n
          a(i,j)=c*a(i,j)
        enddo
      enddo
      return
      end

      subroutine ludcmp(np,n,a,indx,d,det,ierror)
c
c LU decomposition of the matrix a(np,np)
c n   = dimension
c det = determinant
c
      implicit double precision(a-h,o-z)
      parameter (nmax=2000,tiny=1.0d-30)
      dimension a(np,np),indx(n),rvv(nmax)
      d=1.d0
      det=1.d0
      ierror=0
      do 12 i=1,n
        raamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.raamax) raamax=abs(a(i,j))
11      continue
        if (raamax.eq.0.d0)then
        write (*,*) 'singular matrix in ludcmp'
        ierror=1
        return
        endif
        rvv(i)=1.d0/raamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        raamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          rdum=rvv(i)*abs(sum)
          if (rdum.ge.raamax) then
            imax=i
            raamax=rdum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          rvv(imax)=rvv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(abs(a(j,j)).eq.0.d0) a(j,j)=tiny
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(abs(a(n,n)).eq.0.d0) a(n,n)=tiny
      do i=1,n
        det=det*a(i,i)
      enddo
      det=det*d
      return
      end

      subroutine lubksb(np,n,a,indx,b)
c
c solve a*x=b, x--->b
c
      implicit double precision(a-h,o-z)
      dimension a(np,np),indx(n),b(n)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
      return
      end

