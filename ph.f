      call tpb
      stop 
      end


      subroutine tpb

c GLOSSARY:
c______________________________________________________________________
c
c
c______________________________________________________________________


c      implicit double precision(a-h,o-z)
      implicit none
      integer i,j,nl,i1,i2,i3,i4,i5,i6,i7,i8,ns,nums,np
      integer :: ispb_n(8),ispb_s(8),mps(8),b4bd(4,6),b2bd(2,28)

      nl=4
      ns=2*nl
      nums=0
      np  =0

      i=0
         do j=1,nl
            i=i+1
            ispb_n(i)=j-1
            ispb_s(i)=1
            i=i+1
            ispb_n(i)=j-1
            ispb_s(i)=-1
         enddo

c      do i=1,8
c         write(*,*)ispb_n(i),ispb_s(i)
c      enddo

      !do i=1,8
         mps(:)=0
      !enddo

      do i1=1,ns
         mps(ns)=0
         do i2=1,ns
            if(i2.eq.i1) mps(i2)=1
         enddo
         if(i1.gt.1) mps(i1-1)=0
c*** WRITE BLOCK
c         write(*,*)(mps(i),i=1,8)
c         nums=nums+1
c*** WRITE BLOCK
c TWO PARTICLE STATES
         do i3=i1+1,ns
            do i4=1,ns
               if(i4.eq.i3) mps(i4)=1
            enddo
            do i=i1+1,i3-1
               mps(i)=0
            enddo
c*** WRITE BLOCK
c         write(*,*)(mps(i),i=1,8)
c         nums=nums+1
c*** WRITE BLOCK
c THREE PARTICLE STATE   
            do i5=i3+1,ns
               do i6=1,ns
                  if(i6.eq.i5) mps(i6)=1
               enddo
               do i=i3+1,i5-1
                  mps(i)=0
               enddo
c*** WRITE BLOCK
c               do i=1,8
c                  np=np+mps(i)
c               enddo
c               write(*,*)(mps(i),i=1,8),np
c               np=0
c               nums=nums+1
c*** WRITE BLOCK
c               mps(ns)=0
c FOUR PARTICLE STATE
               do i7=i5+1,ns
                  do i8=1,ns
                     if(i7.eq.i8) mps(i8)=1
                  enddo
                  do i=i5+1,i7-1
                     mps(i)=0
                  enddo
c*** WRITE BLOCK
                  do i=1,8
                     np=np+mps(i)
                  enddo
                  write(*,*)(mps(i),i=1,8),np
                  np=0
                  nums=nums+1
c*** WRITE BLOCK
               mps(ns)=0
               enddo
            enddo
         enddo
      write(*,*)''
      enddo
      write(*,*)nums

c SMART WAY TO DO IT
      j=0
      b4bd(:,:)=0
      do i1 = 1 , ns
       do i2 = i1+1 , ns
        do i3 = i2+1 , ns
         do i4 = i3+1 , ns
          mps(:) = 0
          mps(i1) = 1
          mps(i2) = 1
          mps(i3) = 1
          mps(i4) = 1
          if(ispb_s(i1).ne.-1*ispb_s(i2)) goto 10
          if(ispb_n(i1).ne.ispb_n(i2)) goto 10
          if(ispb_s(i3).ne.-1*ispb_s(i4)) goto 10
          if(ispb_n(i3).ne.ispb_n(i4)) goto 10
          write(*,*)(mps(i),i=1,8)
          j=j+1
          b4bd(1,j)=i1
          b4bd(2,j)=i2
          b4bd(3,j)=i3
          b4bd(4,j)=i4
   10 continue
      enddo
      enddo
      enddo
      enddo
      write(*,*)

      do i=1,4
         write(*,*)(b4bd(i,j),j=1,6)
      enddo
      write(*,*)'?!?!'

c TWO-BODY CONFIGURATIONS

      return
      end
      
      

      
      
      
      


