module pbc

   use const, only : PREC, L_INT

   implicit none

   logical, save :: flg_pbc

   ! Current box size
   real(PREC), save :: pbc_box(3)
   real(PREC), save :: pbc_box_half(3)

   ! Store values from the input file
   real(PREC), save :: pbc_box_input(3)
   logical, save :: flg_pbc_ignore_rst

   ! Resizing  box
   logical, save :: flg_pbc_resize
   logical, save :: flg_pbc_resize_target
   integer(L_INT), save :: pbc_resize_step
   integer(L_INT), save :: pbc_resize_count(3)
   real(PREC), save :: pbc_box_original(3)
   real(PREC), save :: pbc_resize_change(3)
   real(PREC), save :: pbc_resize_target(3)

contains

   subroutine init_pbc(restarted)
      use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, INT64
      use const_idx, only : RSTBLK
      logical, intent(in) :: restarted

      integer :: i
      integer :: rst_status

      print '(a)', 'Setting up the Periodic Boundary Conditon'
      print '(a, 3(1x,g15.8))', 'PBC box size in the input file:  ', pbc_box_input(:)
      pbc_box(:) = pbc_box_input(:)

      if (restarted .and. .not. flg_pbc_ignore_rst) then
         call read_rst(RSTBLK%PBC, rst_status)
         ! If the restart file conatins PBC block, pbc_box will be overwritten.
         if (rst_status == 0) then
            print '(a, 3(1x,g15.8))', 'PBC box size in the restart file:', pbc_box(:)
         else
            print '(a)', 'PBC record was not found in the restart file.'
         endif
      endif

      print '(a, 3(1x,g15.8))', 'Set the PBC box size to          ', pbc_box(:)
      pbc_box_half(:) = 0.5_PREC * pbc_box(:)

      ! Initialise values for resize
      pbc_resize_count(1:3) = 0
      pbc_box_original(:) = pbc_box(:)

      if (flg_pbc_resize) then
         print '(a, i12, a)',        'Resizing scheme enabled: For every ', pbc_resize_step, ' steps,'
         print '(a, 3(1x,g15.8),a)', 'the PBC box will be resized by   ', pbc_resize_change(:)
         if (flg_pbc_resize_target) then
            do i = 1, 3
               if (pbc_resize_change(i) > 0.0_PREC) then
                  if (pbc_resize_target(i) < pbc_box(i)) then
                     print *, 'Error: inconsistent setup for PBC resize, target < initial size while the change > 0.'
                     call sis_abort()
                  endif
               else
                  if (pbc_resize_target(i) > pbc_box(i)) then
                     print *, 'Error: inconsistent setup for PBC resize, target > initial size while the change < 0.'
                     call sis_abort()
                  endif
               endif
            enddo
            print '(a, 3(1x,g15.8),a)', 'until reaching the target size,  ', pbc_resize_target(:)
         endif
      endif

      print *

   endsubroutine init_pbc

   subroutine pbc_resize()
      use const_phys, only : SMALL_VALUE 
      integer :: i
      real(PREC) :: x
      logical :: reach_target(3)

      reach_target(:) = .False.

      do i = 1, 3
         ! Do nothing if change is zero
         if ((abs(pbc_resize_change(i)) < SMALL_VALUE)) then
            reach_target(i) = .True.
            cycle
         endif

         ! Do nothing if the current size is already the target
         if (flg_pbc_resize_target) then
            if (abs(pbc_box(i) - pbc_resize_target(i)) < SMALL_VALUE) then
               reach_target(i) = .True.
               cycle
            endif
         endif

         pbc_resize_count(i) = pbc_resize_count(i) + 1
         x = pbc_box_original(i) + pbc_resize_change(i) * pbc_resize_count(i)

         ! Compress
         if (pbc_resize_change(i) < 0.0_PREC) then
            if (flg_pbc_resize_target) then
               if (x - pbc_resize_target(i) < SMALL_VALUE) then
                  pbc_box(i) = pbc_resize_target(i)
                  reach_target(i) = .True.
               else
                  pbc_box(i) = x
               endif
            else
               pbc_box(i) = x
            endif

         else ! Expand
            if (flg_pbc_resize_target) then
               if (x - pbc_resize_target(i) > SMALL_VALUE) then
                  pbc_box(i) = pbc_resize_target(i)
                  reach_target(i) = .True.
               else
                  pbc_box(i) = x
               endif
            else
               pbc_box(i) = x
            endif
         endif
      enddo

      pbc_box_half(:) = 0.5_PREC * pbc_box(:)

      if (all(reach_target)) then
         flg_pbc_resize = .False.
      endif

   endsubroutine pbc_resize

   pure function pbc_vec_d(v1, v2) result (new_vec)

      use const, only : PREC

      real(PREC) :: new_vec(3)
      real(PREC), intent(in) :: v1(3), v2(3)

      integer :: i

      new_vec(:) = v1(:) - v2(:)

      if (.not. flg_pbc) then
         return
      endif

      do i = 1, 3
         if(new_vec(i) > pbc_box_half(i)) then
            new_vec(i) = new_vec(i) - pbc_box(i)

         else if(new_vec(i) < -pbc_box_half(i)) then
            new_vec(i) = new_vec(i) + pbc_box(i)

         end if
      end do

   end function pbc_vec_d

   pure function pbc_vec(v) result (new_vec)
      
      use const, only : PREC

      real(PREC) :: new_vec(3)
      real(PREC), intent(in) :: v(3)

      integer :: i

      if (.not. flg_pbc) then
         new_vec(1:3) = v(1:3)
         return
      endif

      do i = 1, 3
         if(v(i) > pbc_box_half(i)) then
            new_vec(i) = v(i) - pbc_box(i)
   
         else if(v(i) < -pbc_box_half(i)) then
            new_vec(i) = v(i) + pbc_box(i)
   
         else
            new_vec(i) = v(i)
   
         end if
      end do

   end function pbc_vec

   subroutine pbc_wrap(irep)
      
      use const
      use var_top, only : nmp
      use var_state, only : xyz

      implicit none
  
      integer, intent(in) :: irep

      integer :: imp, i
      real(PREC) :: box_max(3), box_min(3)
      real(PREC) :: x
     
      !box_max(1:3) =  0.5*pbc_box(1:3)
      !box_min(1:3) = -0.5*pbc_box(1:3)
      box_max(1:3) = pbc_box(1:3)
      box_min(1:3) = 0.0_PREC
  
      do imp = 1, nmp
         do i = 1, 3
            x = xyz(i, imp, irep)
            if(x > box_max(i)) then
               xyz(i, imp, irep) = x - pbc_box(i) * (int((x - box_max(i))/pbc_box(i)) + 1)
            else if(x < box_min(i)) then
               xyz(i, imp, irep) = x + pbc_box(i) * (int((box_min(i) - x)/pbc_box(i)) + 1)
            endif
         enddo
      enddo

   endsubroutine pbc_wrap


!   subroutine pbc_wrap_chain(irep)
!
!   TO be implemented
!
!   endsubroutine pbc_wrap_chain


endmodule pbc
