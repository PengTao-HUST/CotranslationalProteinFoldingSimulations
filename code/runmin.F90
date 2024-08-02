#include "copyright.i"

!*******************************************************************************
!
! Module: runmin_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module runmin_mod

  implicit none

  private       grdmax, printe

! TIMINGS NOTE - Time not allocated to force determination or mpi communications
!                gets lumped as "other" time rather than being called out as
!                a minimization time; in MD, this would be called runmd_time...

contains

!*******************************************************************************
!
! Subroutine:  runmin_master
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine runmin_master(atm_cnt, crd, frc, igraph)

  use prmtop_dat_mod
  use constraints_mod
  use pmemd_lib_mod
  use constraints_mod
  use extra_pnts_nb14_mod
  use file_io_mod
  use file_io_dat_mod
  use gb_force_mod
  use img_mod
#ifdef MPI
  use loadbal_mod
#endif
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use parallel_mod
  use gb_parallel_mod
  use pme_force_mod
  use shake_mod
  use state_info_mod
  use timers_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3 * atm_cnt)
  double precision      :: frc(3 * atm_cnt)
  character(4)          :: igraph(atm_cnt)

! Local variables:

  logical               :: belly
  logical               :: not_first_loop
  logical               :: new_list
  
  double precision      :: betax
  double precision      :: crits
  double precision      :: ddspln
  double precision      :: dfpr
  double precision      :: dfpred
  double precision      :: dxst
  double precision      :: dxsth
  double precision      :: dxstm
  double precision      :: f
  double precision      :: fch
  double precision      :: fdmax
  double precision      :: finit
  double precision      :: fmin
  double precision      :: fnq
  double precision      :: fold
  double precision      :: gama
  double precision      :: gamden
  double precision      :: ginit
  double precision      :: gmin
  double precision      :: gnew
  double precision      :: gspln
  double precision      :: gsqrd
  integer               :: i
  integer               :: iatmax
  integer               :: ier
  integer               :: ilmnfl
  integer               :: iretry
  integer               :: iterc
  integer               :: iterfm
  integer               :: iterrs
  integer               :: itmout
  integer               :: kstcyc
  character(4)          :: labmax
  integer               :: linmin
  integer               :: maxlin
  integer               :: mstcyc
  integer               :: mxfcon
  logical               :: newstr
  integer               :: atm_cnt3
  integer               :: ncalls
  integer               :: nct
  integer               :: ndfp
  integer               :: nfbeg
  integer               :: nfopt
  integer               :: nstcyc
  integer               :: ntnb
  double precision      :: rms
  double precision      :: sbound
  logical               :: steep
  double precision      :: step
  double precision      :: stepch
  double precision      :: stmin
  double precision      :: sum
  double precision      :: swork
  double precision      :: work
  type(gb_pot_ene_rec)  :: gb_pot_ene
  type(pme_pot_ene_rec) :: pme_pot_ene
  double precision      :: si(si_cnt)   ! system state info array.
  double precision      :: virial(3)    ! Dummy space; returned value ignored.
  double precision      :: ekcmt(3)     ! Dummy space; returned value ignored.
  double precision      :: pme_err_est  ! Dummy space; returned value ignored.
! TI variables:
  integer               :: ti_nct
  integer               :: ib, jb
  
! These are the work arrays, dynamically allocated on the stack:

  double precision      :: w(3 * atm_cnt)
  double precision      :: w_xopt(3 * atm_cnt)
  double precision      :: w_gopt(3 * atm_cnt)
  double precision      :: w_ginit(3 * atm_cnt)
  double precision      :: w_rsdg(3 * atm_cnt)
  double precision      :: w_rsdx(3 * atm_cnt)

  data maxlin, mxfcon, kstcyc / 10, 4 , 4 /
  data dxstm, crits, dfpred /1.0d-05, 1.0d-06, 1.d0/

  call zero_time()

  si(:) = 0.d0

! Evaluate some constants:

  fmin = 0.0d0
  atm_cnt3 = atm_cnt * 3
  belly = ibelly .gt. 0
  ier = 0
  nct = 0
  if (ntc .eq. 2 .or. ntc .eq. 4) nct = nbonh
  if (ntc .eq. 3) nct = nbonh + nbona
  ndfp = atm_cnt3 - nct
  if (belly) ndfp = 3 * belly_atm_cnt - nct
  ntnb = 1 

  if (ti_mode .ne. 0) then
    ! Correct dof in this case
    ti_nct = 0
    ! TI is not compatible with ntc .eq. 3
    if (ntc .eq. 2) then
      do i = 1, nbonh
        ib = gbl_bond(i)%atm_i
        jb = gbl_bond(i)%atm_j
        if (ti_lst(2,ib)+ti_lst(2,jb) .ne. 0 &
            .and. ti_sc_lst(ib)+ti_sc_lst(jb) .eq. 0) cycle
        ti_nct = ti_nct + 1
      end do
    end if
    ndfp = (atm_cnt - ti_latm_cnt(1)) * 3 - ti_nct
  end if

  fnq = sqrt(dble(ndfp))
  rms = 0.0d0
  steep = .false.
  newstr = .false.
  nstcyc = 0
  mstcyc = kstcyc
  if (ntmin .eq. 2) mstcyc = maxcyc
  if (ntmin .eq. 1) mstcyc = ncyc
  if (ntmin .gt. 0) steep = .true.
  fold = 0.0d0
  dxst = dx0
  linmin = 0
  itmout = 0
  ilmnfl = 0

! Set some parameters to begin the calculation:

  iterc = 0
  ncalls = 0
  iterfm = iterc

! Let the initial search direction be minus the gradient vector. iterrs gives
! the iteration number of the most recent restart, but is set to zero when
! steepest descent direction is used:

! (Here is the beginning of a big loop:)

  not_first_loop = .false.
#ifdef MPI
  notdone = 1
#endif

   20 continue

! Gather the submolecules into the box:

  ncalls = ncalls + 1
  if (ncalls .eq. ncalls/nsnb*nsnb) ntnb = 1
  if (ntnb .eq. 1 .and. ncalls .gt. 1) steep = .true.

! Calculate the force and energy:

  call update_time(other_time)

! Apply shake to constrain bonds if necessary:

  if (ntc .ne. 1) then
    frc(:) = crd(:)
    call shake(frc, crd)
    call shake_fastwater(frc, crd)
    call update_time(shake_time)
  end if

! If using extra points, update them here, after all coordinate movement in
! shake or in the previous cycle of the minimization code...

  if (numextra .gt. 0 ) then
    if (frameon .ne. 0 .and. gbl_frame_cnt .ne. 0) &
      call all_local_to_global(crd, ep_frames, ep_lcl_crd, gbl_frame_cnt)
  end if


#ifdef MPI
  call mpi_bcast(notdone, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

  ! For minimization, everyone gets a full set of crds from the master;
  ! minimization is unfortunately not well parallelized.

  call mpi_bcast(crd, 3 * atm_cnt, mpi_double_precision, 0, & 
                 pmemd_comm, err_code_mpi)
  call update_time(fcve_dist_time)
#endif

  if (using_pme_potential) then

    if (not_first_loop) then

      ! Now do a skin check to see if we will have to rebuild the pairlist.

      call check_all_atom_movement(atm_cnt, crd, gbl_atm_saved_crd, skinnb, &
                                   ntp, new_list)

#ifdef MPI
      call check_new_list_limit(new_list)
#endif
      call update_time(nonbond_time)
    else

      new_list = .true.
      not_first_loop = .true.

    end if

#ifdef CUDA
    call gpu_upload_crd(crd)
#endif
#ifdef MPI
    call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
                   gbl_my_atm_lst, new_list, .true., .false., pme_pot_ene, &
                   ncalls, virial, ekcmt, pme_err_est)
#else
    call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
                   new_list, .true., .false., pme_pot_ene, &
                   ncalls, virial, ekcmt, pme_err_est)
#endif
#ifdef CUDA
    call gpu_download_frc(frc)
#endif

    ! Store energy terms in state info array for printout.

    si(si_pot_ene) = pme_pot_ene%total
    si(si_vdw_ene) = pme_pot_ene%vdw_tot
    si(si_elect_ene) = pme_pot_ene%elec_tot
    si(si_hbond_ene) = pme_pot_ene%hbond
    si(si_bond_ene) = pme_pot_ene%bond
    si(si_angle_ene) = pme_pot_ene%angle
    si(si_dihedral_ene) = pme_pot_ene%dihedral
    si(si_vdw_14_ene) = pme_pot_ene%vdw_14
    si(si_elect_14_ene) = pme_pot_ene%elec_14
    si(si_restraint_ene) = pme_pot_ene%restraint
    si(si_angle_ub_ene) = pme_pot_ene%angle_ub
    si(si_dihedral_imp_ene) = pme_pot_ene%imp
    si(si_cmap_ene) = pme_pot_ene%cmap

  else if (using_gb_potential) then
#ifdef CUDA
    call gpu_upload_crd(crd)
#endif
    !always calculate pot_enes when doing minimization - may be overkill?
    call gb_force(atm_cnt, crd, frc, gb_pot_ene, ncalls, .true.)
#ifdef CUDA
    call gpu_download_frc(frc)
#endif

    si(si_pot_ene) = gb_pot_ene%total
    si(si_vdw_ene) = gb_pot_ene%vdw_tot
    si(si_elect_ene) = gb_pot_ene%elec_tot
    si(si_hbond_ene) = gb_pot_ene%gb                ! temporary hack
    si(si_surf_ene) = gb_pot_ene%surf
    si(si_bond_ene) = gb_pot_ene%bond
    si(si_angle_ene) = gb_pot_ene%angle
    si(si_dihedral_ene) = gb_pot_ene%dihedral
    si(si_vdw_14_ene) = gb_pot_ene%vdw_14
    si(si_elect_14_ene) = gb_pot_ene%elec_14
    si(si_restraint_ene) = gb_pot_ene%restraint
    si(si_pme_err_est) = 0.d0
    si(si_angle_ub_ene) = gb_pot_ene%angle_ub
    si(si_dihedral_imp_ene) = gb_pot_ene%imp
    si(si_cmap_ene) = gb_pot_ene%cmap

  end if

#ifdef MPI
  if (using_pme_potential) then
    call mpi_gathervec(atm_cnt, frc)
  else if (using_gb_potential) then
    call gb_mpi_gathervec(atm_cnt, frc)
  end if
  call update_time(fcve_dist_time)
#endif

  if (ti_mode .ne. 0) then
    if (ti_latm_cnt(1) .gt. 0) then
#ifdef GTI
    call gti_sync_vector(0, 0)! mixing forces
    call gpu_download_frc(frc)
#else      
      call ti_exchange_vec(atm_cnt, frc, .true.)
#endif      
    end if
  end if

  f = si(si_pot_ene)

  if (ti_mode .ne. 0) then
    si(si_dvdl) = ti_ene(1, si_dvdl) 
    f = f + ti_ene(1, si_pot_ene) + ti_ene(2, si_pot_ene)
  end if

  ntnb = 0
  sum = dot_product(frc, frc)
  rms = sqrt(sum)/fnq

! Print the intermediate results:

  if (mod(ncalls, ntpr) .eq. 0  .or. ncalls .eq. 1) then
    if (master) then
       call grdmax(atm_cnt3, frc, iatmax, fdmax)
       iatmax = (iatmax - 1)/3 + 1
       labmax = igraph(iatmax)
       if (ti_mode .eq. 0) then
         call printe(ncalls, rms, fdmax, si, iatmax, labmax)
         if (nmropt .ne. 0) call nmrptx(6)
       else
         do ti_region = 1, 2
           write(mdout, 601) ti_region       
           call printe(ncalls, rms, fdmax, si, iatmax, labmax)
           if (nmropt .ne. 0) call nmrptx(6)
           call ti_print_ene(ti_region, ti_ene(ti_region,:), &
                             ti_ene_aug(ti_region,:))
         end do 
       end if
    end if
  end if

! Do some steepest steps before entering the conjugate gradient method:

  if (steep) then
    nstcyc = nstcyc + 1
    if (nstcyc .le. mstcyc) then

      if (dxst .le. crits) dxst = dxstm
      dxst = dxst/2.0d0
      if (f .lt. fold) dxst = dxst*2.4d0
      dxsth = dxst/sqrt(sum)
      if (nstcyc .le. 1 .or. f .le. fmin) then
        fmin = f
        nfopt = ncalls
        w_xopt(:) = crd(:)
        w_gopt(:) = - frc(:)
      end if

! Check for convergence:

      if (rms .le. drms) goto 300

      if (ncalls .ge. maxcyc) then
        ier = 131
        goto 290
      end if

      fold = f
      crd(:) = crd(:) + dxsth * frc(:)

      goto 20

    else

! (arrive here when finished with this set of steepest descent cycles)

      steep = .false.
      newstr = .true.
      nstcyc = 0
      mstcyc = kstcyc
    end if
  end if

! Start of conjugate gradient steps:

  frc(:) = - frc(:)

  if (.not. newstr .and. ncalls .ge. 2) goto 82

70 continue

  w(:) = -frc(:)

  iterrs = 0
  if (newstr) iterc = 0

  if (iterc .gt. 0) goto 140

82 continue

  gnew = dot_product(w, frc)

  if (newstr .or. ncalls .eq. 1) goto 100

  fch = f - fmin

! Store the values of crd, f and g, if they are the best that have been
! calculated so far. Test for convergence:

  if (fch) 100, 90, 130

90 continue

  if (gnew/gmin.lt. - 1.0d0) goto 120

100 continue

  fmin = f
  gsqrd = sum
  nfopt = ncalls

  w_xopt(:) = crd(:)
  w_gopt(:) = frc(:)

120 continue

  if (rms .le. drms) goto 300

! Test if the value of maxcyc allows another call of funct:

130 continue

  if (ncalls .ge. maxcyc) then
    ier = 131
    goto 290
  end if

  if (.not.newstr .and. ncalls .gt. 1) goto 180

! This section is executed at the beginning of a conjugate gradient set of
! minimization steps.

! Set dfpr to dx0*gsqrd. dfpr is the reduction in the function value. stmin is
! usually the step-length of the most recent line search that gives the least
! value of f:

! dac change, 10/91:  Return to original idea of trying to go downhill by the
!                     absolute amount, dfpred (which defaults to 1 kcal/mol,
!                     see data statement above).  This can eliminate very bad
!                     initial conjugate gradient steps.

  dfpr = dfpred
  stmin = dfpred/gsqrd

  newstr = .false.

! Begin the main conguate gradient iteration:

140 continue

  iterc = iterc + 1

  finit = f
  ginit = 0.0d0

  w_ginit(:) = frc(:)

  ginit = dot_product(w, frc)
  if (ginit .ge. 0.0d0) goto 260
  gmin = ginit
  sbound = - 1.0d0
  nfbeg = ncalls
  iretry = - 1

  stepch = min(stmin, abs(dfpr/ginit))
  stmin = dxstm

160 continue

  step = stmin + stepch
  dxst = step
  swork = 0.0d0

  crd(:) = w_xopt(:) + stepch*w(:)

  do i = 1, atm_cnt3
    swork = max(swork, abs(crd(i) - w_xopt(i)))
  end do

  if (swork .gt. 0.0d0) goto 20

! "work = swork" may not be needed - wont hurt.  -gls

  work = swork

! Terminate the line search if stepch is effectively zero:

  if (ncalls .gt. nfbeg + 1 .or. abs(gmin/ginit) .gt. 0.2d0) then
    if (master) write(mdout, 370)
    steep = .true.
    linmin = linmin + 1
  end if

  goto 270

180 continue

  work = (fch + fch)/stepch - gnew - gmin
  ddspln = (gnew - gmin)/stepch
  if (ncalls .gt. nfopt) then
    sbound = step
  else
    if (gmin*gnew .le. 0.0d0) sbound = stmin
    stmin = step
    gmin = gnew
    stepch = - stepch
  end if

  if (fch .ne. 0.0d0) ddspln = ddspln + (work + work)/stepch

! Test for convergence of the line search, but force atleast two steps to be
! taken in order not to lose quadratic termination:

  if (gmin .eq. 0.0d0) goto 270
  if (ncalls .le. nfbeg + 1) goto 200
  if (abs(gmin/ginit) .le. 0.2d0) goto 270

! Apply the test that depends on the parameter maxlin:

190 continue

  if (ncalls .lt. nfopt + maxlin) goto 200

! Possible non bonded update. make a restart:

  if (master) write(mdout, 370)
  steep = .true.
  linmin = linmin + 1

  goto 270

200 continue

  stepch = 0.5d0*(sbound-stmin)
  if (sbound.lt. - 0.5d0) stepch = 9.0d0*stmin
  gspln = gmin + stepch*ddspln
  if (gmin*gspln .lt. 0.0d0) stepch = stepch*gmin/(gmin - gspln)

  goto 160

! Calculate the value of betax in the new direction:

210 continue

  sum = dot_product(frc, w_ginit)
  betax = (gsqrd-sum)/(gmin - ginit)

! Test that the new search direction can be made downhill.  If not then try to
! improve the accuracy of the line search:

  if (abs(betax*gmin) .le. 0.2d0*gsqrd) goto 220
  iretry = iretry + 1
  if (iretry .le. 0) goto 190

220 continue

  if (f .lt. finit) iterfm = iterc
  if (iterc .ge. iterfm + mxfcon) then
    if (master) write(mdout, 370)
    steep = .true.
    linmin = linmin + 1
    goto 270
  end if
  dfpr = stmin*ginit

! Branch if a restart procedure is required due to the iteration number or due
! to the scalar product of consecutive gradients:

  if (iretry .gt. 0) goto 70
  if (iterrs .eq. 0) goto 240
  if (iterc - iterrs .ge. atm_cnt3) goto 240
  if (abs(sum) .ge. 0.2d0*gsqrd) goto 240

! Calculate gama in the new search direction. gamden is set by the restart
! procedure:

  gama = dot_product(frc, w_rsdg)
  sum  = dot_product(frc, w_rsdx)
  gama = gama/gamden

! Restart if the new search direction is not sufficiently downhill:

  if (abs(betax*gmin + gama*sum) .ge. 0.2d0*gsqrd) goto 240

! Calculate the new search direction:

  w(:) = - frc(:) + betax*w(:) + gama*w_rsdx(:)

! Cycle back for more conjugate gradient steps:

  goto 140

! Apply the restart procedure:

240 continue

  gamden = gmin - ginit

  w_rsdx(:) = w(:)
  w_rsdg(:) = frc(:) - w_ginit(:)
  w(:) = - frc(:) + betax*w(:)

  iterrs = iterc

  goto 140

! Set ier to indicate that the search direction is uphill:

260 continue

  steep = .true.
  if (master) write(mdout, 370)
  linmin = linmin + 1

! Ensure that f, crd and g are optimal:

270 continue

  if (ncalls .ne. nfopt) then
    f = fmin

    crd(:) = w_xopt(:)
    frc(:) = w_gopt(:)

  end if

  if (linmin .gt. 4) then
    ilmnfl = 1
    goto 310
  end if

  if (steep) goto 20
  if (ier .eq. 0) goto 210

290 continue

  if (master) then
    if (ier .eq. 129) write(mdout, 320)
    if (ier .eq. 130) write(mdout, 330)
    if (ier .eq. 131) write(mdout, 340)
    if (ier .eq. 132) write(mdout, 350)
  end if

  goto 310

300 continue

! Write the final results:

310 continue

  if (master) then
    write(mdout, 380)
    call grdmax(atm_cnt3, frc, iatmax, fdmax)
    iatmax = (iatmax - 1)/3 + 1
    labmax = igraph(iatmax)
    if (ti_mode .eq. 0) then
      call printe(ncalls, rms, fdmax, si, iatmax, labmax)
      if (nmropt .ne. 0) then
        call nmrptx(6)
        call ndvptx(crd, frc, 6)
      end if
    else
      do ti_region = 1, 2
        write(mdout, 601) ti_region       
        call printe(ncalls, rms, fdmax, si, iatmax, labmax)
        if (nmropt .ne. 0) then
          call nmrptx(6)
          call ndvptx(crd, frc, 6)
        end if
        call ti_print_ene(ti_region, ti_ene(ti_region,:),&
                          ti_ene_aug(ti_region,:))
      end do    
    end if

    if (itmout .eq. 1) write(mdout, 360)
    if (ilmnfl .eq. 1) write(mdout, 390)
  end if

  return

320 format('  LINE SEARCH ABANDONED ... PROBLEM WITH G')
330 format('  SEARCH DIRECTION IS UPHILL ')
340 format(/,/'  Maximum number of minimization cycles reached.')
350 format('  THE VALUE OF F COULD NOT BE REDUCED')
360 format(/4x, 'WARNING ... TIME LIMIT EXCEEDED ... TO BE RESTARTED')
370 format(/4x, ' .... RESTARTED DUE TO LINMIN FAILURE ...')
380 format(/ /20x, 'FINAL RESULTS', /)
390 format(/5x, '***** REPEATED LINMIN FAILURE *****')
601 format(/,'| TI region ', i2, /)

end subroutine runmin_master

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  runmin_slave
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine runmin_slave(atm_cnt, crd, frc, img_atm_map)

  use gb_force_mod
  use pmemd_lib_mod
  use img_mod
#ifdef MPI
  use loadbal_mod
#endif
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use parallel_dat_mod
  use parallel_mod
  use gb_parallel_mod
  use pme_force_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3 * atm_cnt)
  double precision              :: frc(3 * atm_cnt)
  integer                       :: img_atm_map(atm_cnt)

! Local variables:

  integer                       :: ncalls
  logical                       :: new_list
  logical                       :: not_first_loop
  type(gb_pot_ene_rec)          :: gb_pot_ene
  type(pme_pot_ene_rec)         :: pme_pot_ene

  double precision      :: virial(3)    ! Dummy space; returned value ignored.
  double precision      :: ekcmt(3)     ! Dummy space; returned value ignored.
  double precision      :: pme_err_est  ! Dummy space; returned value ignored.

  call zero_time()

  not_first_loop = .false.

  ncalls = 0

  do

    call mpi_bcast(notdone, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
    call update_time(fcve_dist_time)

    if (notdone .ne. 1) return

    ! For minimization, everyone gets a full set of crds from the master;
    ! minimization is unfortunately not well parallelized.

    call mpi_bcast(crd, 3 * atm_cnt, mpi_double_precision, 0, & 
                   pmemd_comm, err_code_mpi)
    call update_time(fcve_dist_time)

    if (using_pme_potential) then

      if (not_first_loop) then

        ! Do a skin check to see if we will have to rebuild the pairlist.

        call check_all_atom_movement(atm_cnt, crd, gbl_atm_saved_crd, skinnb, &
                                     ntp, new_list)
#ifdef MPI
        call check_new_list_limit(new_list)
#endif
        call update_time(nonbond_time)
      else

        new_list = .true.
        not_first_loop = .true.

      end if

#ifdef MPI
      call pme_force(atm_cnt, crd, frc, img_atm_map, gbl_atm_img_map, &
                     gbl_my_atm_lst, new_list, .true., .false., pme_pot_ene, &
                     ncalls, virial, ekcmt, pme_err_est)
#else
      call pme_force(atm_cnt, crd, frc, img_atm_map, gbl_atm_img_map, &
                       new_list, .true., .false., pme_pot_ene, &
                       ncalls, virial, ekcmt, pme_err_est)
#endif /* MPI */
    else if (using_gb_potential) then

      ncalls = ncalls + 1
      !always calculate pot_enes when doing minimization - may be overkill -
      !needed for mpi slaves?
      call gb_force(atm_cnt, crd, frc, gb_pot_ene, ncalls, .true.)

    end if

    ! Potential energy info not used in slaves...

    if (using_pme_potential) then
      call mpi_gathervec(atm_cnt, frc)
    else if (using_gb_potential) then
      call gb_mpi_gathervec(atm_cnt, frc)
    end if
    call update_time(fcve_dist_time)

  end do

  return

end subroutine runmin_slave
#endif

!*******************************************************************************
!
! Internal Subroutine:  grdmax
!
! Description: <TBS>
!              
! Rev A mods:  Converted this routine from function to subrt.
!              Added iatmax = atom number of max gradient to args.
!
!*******************************************************************************

subroutine grdmax(n, g, iatmax, fdmax)

  implicit none

! Formal arguments:

  integer           n
  double precision  g(*)
  integer           iatmax
  double precision  fdmax

! Local variables:

  double precision  dum
  integer           i
  double precision  gi

  dum = 0.0d0
  iatmax = 1
  do i = 1, n
    gi = abs(g(i))
    if (gi .gt. dum) then
      dum = gi
      iatmax = i
    end if
  end do
  fdmax = dum

  return

end subroutine grdmax

!*******************************************************************************
!
! Internal Subroutine:  printe
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine printe(nstep, dff, fdmax, si, iatmax, labmax)

  use file_io_mod
  use file_io_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod
  use pme_force_mod
  use runfiles_mod
  use state_info_mod
  use charmm_mod, only : charmm_active
  use ti_mod
  implicit none

! Formal arguments:

  integer                       :: nstep
  double precision              :: dff
  double precision              :: fdmax
  double precision              :: si(si_cnt)
  integer                       :: iatmax
  character(4)                  :: labmax

! Local variables:

  logical, save         :: first_6call = .true.
  integer, save         :: next_6flush_sec
  logical, save         :: first_7call = .true.
  integer, save         :: next_7flush_sec
  integer               :: current_sec
  integer               :: current_usec           ! Dummy, not used.

! Partition energy terms:

  write(mdout, 9018)
  write(mdout, 9028) nstep, si(si_pot_ene), dff, fdmax, labmax, iatmax
  write(mdout, 9038) si(si_bond_ene), si(si_angle_ene), si(si_dihedral_ene)
  if (charmm_active) then
    write(mdout, 9039) si(si_angle_ub_ene), si(si_dihedral_imp_ene), &
                       si(si_cmap_ene)
  endif

  if (using_pme_potential) then
    write(mdout, 9048) si(si_vdw_ene), si(si_elect_ene), &
                       si(si_hbond_ene)
  else
    write(mdout, 9049) si(si_vdw_ene), si(si_elect_ene), &
                       si(si_hbond_ene)
  end if

  write(mdout, 9059) si(si_vdw_14_ene), si(si_elect_14_ene), &
                     si(si_restraint_ene)

  if (gbsa .eq. 1) &
    write(mdout, 9050) si(si_surf_ene)

  if (si(si_restraint_ene) .ne. 0.0) &
    write(mdout, 9078) si(si_pot_ene) - si(si_restraint_ene)
  if (ti_mode .ne. 0) write(mdout, 9100)  si(si_dvdl)
! Check if we need to force a flush of mdout. Barring changes in the unix
! system call, the clock is going to wrap in 2038, and flushing won't be
! strictly correct for a flush_interval...

  call get_wall_time(current_sec, current_usec)

  if (first_6call) then
    first_6call = .false.
    next_6flush_sec = current_sec + mdout_flush_interval
    close(mdout)
    open(unit=mdout, file=mdout_name, status='OLD', position='APPEND')
  else
    if (current_sec .ge. next_6flush_sec) then
      next_6flush_sec = current_sec + mdout_flush_interval
      close(mdout)
      open(unit=mdout, file=mdout_name, status='OLD', position='APPEND')
    end if
  end if

! Flush i/o buffer:

! Flushing actually does not work particularly reliably for a number of
! machines and compilers, and in the more benign cases simply fails, but in
! the more malign cases can actually corrupt the stack (due to a compiler-
! dependent flush() call interface change).  We therefore no longer do
! flushes of anything in PMEMD; if it needs to go out, we close it and reopen
! it.

! Output the mdinfo file if requested, and if the flush interval has elapsed.

  if (first_7call) then
    first_7call = .false.
    next_7flush_sec = current_sec + mdinfo_flush_interval
  else
    if (current_sec .lt. next_7flush_sec) return
    next_7flush_sec = current_sec + mdinfo_flush_interval
  end if

  call amopen(mdinfo, mdinfo_name, 'U', 'F', 'W')

  write(mdinfo, 9018)
  write(mdinfo, 9028) nstep, si(si_pot_ene), dff, fdmax, labmax, iatmax
  write(mdinfo, 9038) si(si_bond_ene), si(si_angle_ene), si(si_dihedral_ene)
  if (charmm_active) then
    write(mdinfo, 9039) si(si_angle_ub_ene), si(si_dihedral_imp_ene), &
                        si(si_cmap_ene)
  endif

  if (using_pme_potential) then
    write(mdinfo, 9048) si(si_vdw_ene), si(si_elect_ene), &
                        si(si_hbond_ene)
  else
    write(mdinfo, 9049) si(si_vdw_ene), si(si_elect_ene), &
                        si(si_hbond_ene)
    if (gbsa .eq. 1) &
       write(mdinfo, 9050) si(si_surf_ene)
  end if

  write(mdinfo, 9059) si(si_vdw_14_ene), si(si_elect_14_ene), &
                     si(si_restraint_ene)

  if (si(si_restraint_ene) .ne. 0.0) &
    write(mdinfo, 9078) si(si_pot_ene) - si(si_restraint_ene)

  if (nmropt .ne. 0) call nmrptx(mdinfo)
  if (ti_mode .ne. 0) write(mdinfo, 9100)  si(si_dvdl)
  close(mdinfo)

  return

9018 format(/ /, 3x, 'NSTEP', 7x, 'ENERGY', 10x, 'RMS', 12x, 'GMAX', 9x, &
            'NAME', 4x, 'NUMBER')
9028 format(1x, i6, 2x, 3(2x, 1pe13.4), 5x, a4, 2x, i7, /)
9038 format(1x, 'BOND    = ', f13.4, 2x, 'ANGLE   = ', f13.4, 2x, &
            'DIHED      = ', f13.4)
9039 format(1x, 'UB      = ', f13.4, 2x, 'IMP     = ', f13.4, 2x, &
            'CMAP       = ', f13.4)
9048 format(1x, 'VDWAALS = ', f13.4, 2x, 'EEL     = ', f13.4, 2x, &
            'HBOND      = ', f13.4)
9049 format(1x, 'VDWAALS = ', f13.4, 2x, 'EEL     = ', f13.4, 2x, &
            'EGB        = ', f13.4)
9050 format (1x,'ESURF   = ',f13.4)
9059 format(1x, '1-4 VDW = ', f13.4, 2x, '1-4 EEL = ', f13.4, 2x, &
            'RESTRAINT  = ', f13.4)
9078 format(1x, 'EAMBER  = ', f13.4)
9100 format (1x,'DV/DL  = ',f14.4)
end subroutine printe

end module runmin_mod
