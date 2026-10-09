!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
!*******************************************************************
!
!  MODULE MOD_TVV
!
!  Time-varying crest vertical element walls (TVV): overflow gates on
!  IBTYPE=64 boundaries whose crest elevation z_c(t) follows a time
!  series given in a crest table. A gate is the set of VEW node pairs
!  whose fort.142 (TVW file) lines name the same crest table
!  (VaryType=4).
!
!  VaryType=4 lines are read by every subdomain from the fulldomain TVW
!  file, so all ranks process the same list of pairs. Each pair is
!  located by a search point (X1, Y1) and a radius; the ranks then
!  agree on the pair through global node numbers.
!
!  2026-10-08: Created by SB
!
!*******************************************************************
module mod_tvv
   use mod_logging, only: allMessage, logMessage, INFO, ERROR, WARNING
   use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
   implicit none
   private

   real(8), parameter :: TVV_NULL = -99999d0
   real(8), parameter :: TVV_DEFAULT_RADIUS = 1d-6 !...m, used when no radius is given
   real(8), parameter :: TVV_TOL = 1d-6 !...m, tolerance of the startup consistency checks

   !...A gate: VEW pairs that follow one crest table
   type :: t_tvv_gate
      character(1024) :: table_file = ''
      integer :: hot = 0 !...1: table times are relative to the hot-start time
      real(8) :: toffset = 0d0 !...s, added to the table times
      integer :: nrec = 0
      real(8), allocatable :: t(:) !...s
      real(8), allocatable :: z(:) !...m, crest elevation
      integer :: npairs = 0
      integer, allocatable :: gtop(:), gbed(:) !...fulldomain node numbers of the wall-top and bed nodes
      integer, allocatable :: ltop(:), lbed(:) !...local node numbers (0 = not in this subdomain)
      integer, allocatable :: etop(:), ebed(:) !...local IBTYPE=64 boundary entries with NBV = top / bed (0 = none)
      real(8), allocatable :: delta(:) !...m, crest minus wall-top elevation
      real(8), allocatable :: zmin(:) !...m, lowest crest: bed elevation + delta
      real(8), allocatable :: zc0(:) !...m, crest in fort.14 (BARINHT)
   end type t_tvv_gate

   logical, public :: tvv_active = .false.
   integer, public :: n_tvv_gates = 0
   type(t_tvv_gate), allocatable, public :: tvv_gates(:)

   public :: tvv_setup, tvv_crest_at

   !...Variables of the TimeVaryingWeir namelist. The group holds every variable
   !   of the namelist read in weir_boundary.F90, so any valid line can be read.
   character(200) :: ScheduleFile
   real(8) :: X1, X2, Y1, Y2
   real(8) :: ZF, ETA_MAX
   real(8) :: TimeStartDay, TimeStartHour, TimeStartMin, TimeStartSec
   real(8) :: TimeEndDay, TimeEndHour, TimeEndMin, TimeEndSec
   real(8) :: FailureDurationDay, FailureDurationHour
   real(8) :: FailureDurationMin, FailureDurationSec
   real(8) :: SearchRadius
   integer :: VaryType, HOT, LOOP, NLOOPS

   namelist /TimeVaryingWeir/ &
      X1, Y1, X2, Y2, &
      VaryType, ZF, &
      ETA_MAX, &
      TimeStartDay, TimeStartHour, TimeStartMin, TimeStartSec, &
      TimeEndDay, TimeEndHour, TimeEndMin, TimeEndSec, &
      FailureDurationDay, FailureDurationHour, &
      FailureDurationMin, FailureDurationSec, &
      ScheduleFile, HOT, LOOP, NLOOPS, SearchRadius

contains

   !-----------------------------------------------------------------------
   !  Read the VaryType=4 lines of the fulldomain TVW file, locate their VEW
   !  pairs, group them into gates and read the crest tables. Called after
   !  WEIR_SETUP at cold and hot start. Does not change any model state.
   !-----------------------------------------------------------------------
   subroutine tvv_setup()
      use SIZES, only: GLOBALDIR
      use GLOBAL, only: TVW_FILE
      use WEIR, only: found_tvw_nml

      integer :: nlines
      real(8), allocatable :: xq(:), yq(:), rq(:)
      integer, allocatable :: line_gate(:), line_no(:)
      character(len=1024), allocatable :: line_table(:)
      integer, allocatable :: gtop(:), gbed(:)
      real(8), allocatable :: delta(:), zbed(:), zc0(:)

      if (.not. found_tvw_nml) return

      call read_tvv_lines(trim(GLOBALDIR)//'/'//trim(TVW_FILE), nlines, xq, yq, rq, &
                          line_table, line_gate, line_no)
      if (nlines == 0) return

      call locate_pairs(nlines, xq, yq, rq, line_no, gtop, gbed, delta, zbed, zc0)
      call build_gates(nlines, line_gate, line_no, gtop, gbed, delta, zbed, zc0)
      call read_tables()
      call log_gates()
      call check_startup()

      tvv_active = n_tvv_gates > 0

   end subroutine tvv_setup

   !-----------------------------------------------------------------------
   !  Crest elevation of gate ig at time timeloc: linear interpolation of
   !  the crest table, held constant outside it. Depends on time only.
   !-----------------------------------------------------------------------
   pure real(8) function tvv_crest_at(ig, timeloc) result(zc)
      integer, intent(in) :: ig
      real(8), intent(in) :: timeloc
      real(8) :: t
      integer :: k

      associate (g => tvv_gates(ig))
         t = timeloc - g%toffset
         if (t <= g%t(1)) then
            zc = g%z(1)
         elseif (t >= g%t(g%nrec)) then
            zc = g%z(g%nrec)
         else
            k = 1
            do while (g%t(k + 1) < t)
               k = k + 1
            end do
            zc = g%z(k) + (g%z(k + 1) - g%z(k))*(t - g%t(k))/(g%t(k + 1) - g%t(k))
         end if
      end associate
   end function tvv_crest_at

   !-----------------------------------------------------------------------
   !  Read the VaryType=4 lines: search point, radius and gate (one gate per
   !  crest table). line_no is the line's position among all entries in the
   !  file, for messages.
   !-----------------------------------------------------------------------
   subroutine read_tvv_lines(fname, nlines, xq, yq, rq, gate_tables, line_gate, line_no)
      use mod_io, only: openFileForRead
      use GLOBAL, only: TVV_SEARCH_RADIUS, IHOT, ITHS, DTDP
      character(*), intent(in) :: fname
      integer, intent(out) :: nlines
      real(8), allocatable, intent(out) :: xq(:), yq(:), rq(:)
      character(len=1024), allocatable, intent(out) :: gate_tables(:)
      integer, allocatable, intent(out) :: line_gate(:), line_no(:)

      integer :: lun, ios, nentries, i, ig, k
      character(2000) :: line
      character(2200) :: nml
      character(512) :: iomsg_
      character(1024) :: msg
      integer, allocatable :: hots(:)
      real(8) :: r

      nlines = 0
      allocate (xq(0), yq(0), rq(0), gate_tables(0), line_gate(0), line_no(0))
      lun = 9142
      call openFileForRead(lun, fname, ios, required=.false.)
      if (ios /= 0) return
      read (lun, *, iostat=ios) nentries
      if (ios /= 0) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="TVV: cannot read the number of entries in "//trim(fname))
      end if

      deallocate (xq, yq, rq, line_gate, line_no)
      allocate (xq(nentries), yq(nentries), rq(nentries), line_gate(nentries), line_no(nentries))
      allocate (hots(nentries))

      do i = 1, nentries
         read (lun, '(A)', iostat=ios) line
         if (ios /= 0) then
            write (msg, '(A,I0,A)') "TVV: cannot read entry ", i, " of "//trim(fname)
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
         end if
         call nullify_nml()
         nml = "&TimeVaryingWeir "//trim(adjustl(line))//" /"
         read (nml, nml=TimeVaryingWeir, iostat=ios, iomsg=iomsg_)
         if (ios /= 0) then
            write (msg, '(A,I0,A)') "TVV: error reading entry ", i, " of "//trim(fname)// &
               ": "//trim(iomsg_)
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
         end if
         if (VaryType /= 4) cycle

         !...Checks of the VaryType=4 line
         write (msg, '(A,I0,A)') "TVV: entry ", i, " (VaryType=4): "
         if (is_null(X1) .or. is_null(Y1)) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg)//" X1= and Y1= are required.")
         end if
         if (.not. (is_null(X2) .and. is_null(Y2))) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg)// &
                           " X2=/Y2= are not used; a pair is located by X1, Y1 and SearchRadius.")
         end if
         if (trim(adjustl(ScheduleFile)) == "NOFILE") then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg)// &
                           " ScheduleFile= (crest table) is required.")
         end if
         if (HOT == -99999) HOT = 0
         if (HOT /= 0 .and. HOT /= 1) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg)//" HOT must be 0 or 1.")
         end if
         if (.not. is_null(SearchRadius)) then
            r = SearchRadius
         elseif (.not. is_null(TVV_SEARCH_RADIUS)) then
            r = TVV_SEARCH_RADIUS
         else
            r = TVV_DEFAULT_RADIUS
         end if
         if (r <= 0d0) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg)//" the search radius must be positive.")
         end if

         nlines = nlines + 1
         call to_model_xy(X1, Y1, xq(nlines), yq(nlines))
         rq(nlines) = r
         line_no(nlines) = i

         !...Gate = crest table; all lines of a gate must agree on HOT
         ig = 0
         do k = 1, size(gate_tables)
            if (trim(gate_tables(k)) == trim(adjustl(ScheduleFile))) then
               ig = k
               exit
            end if
         end do
         if (ig == 0) then
            gate_tables = [character(len=1024) :: gate_tables, trim(adjustl(ScheduleFile))]
            ig = size(gate_tables)
            hots(ig) = HOT
         elseif (hots(ig) /= HOT) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg)// &
                           " all entries naming "//trim(gate_tables(ig))//" must use the same HOT.")
         end if
         line_gate(nlines) = ig
      end do
      close (lun)

      !...Create the gates (pairs are added in build_gates)
      n_tvv_gates = size(gate_tables)
      allocate (tvv_gates(n_tvv_gates))
      do ig = 1, n_tvv_gates
         tvv_gates(ig)%table_file = gate_tables(ig)
         tvv_gates(ig)%hot = hots(ig)
         tvv_gates(ig)%toffset = 0d0
         if (hots(ig) == 1 .and. IHOT /= 0) tvv_gates(ig)%toffset = DTDP*dble(ITHS)
      end do

   end subroutine read_tvv_lines

   !-----------------------------------------------------------------------
   !  Locate the VEW pair of each line. Every rank looks among its own
   !  IBTYPE=64 boundary entries within the search radius; the ranks then
   !  agree on one pair per line through fulldomain node numbers.
   !-----------------------------------------------------------------------
   subroutine locate_pairs(nlines, xq, yq, rq, line_no, gtop, gbed, delta, zbed, zc0)
      use MESH, only: X, Y, DP
      use BOUNDARIES, only: NVEL, LBCODEI, NBV, IBCONN, BARINHT
      integer, intent(in) :: nlines
      real(8), intent(in) :: xq(:), yq(:), rq(:)
      integer, intent(in) :: line_no(:)
      integer, allocatable, intent(out) :: gtop(:), gbed(:)
      real(8), allocatable, intent(out) :: delta(:), zbed(:), zc0(:)

      integer, allocatable :: plo(:), phi(:), mlo(:), mhi(:), seen1(:), seen2(:), flag(:)
      integer :: l, i, n1, n2, g1, g2, lo, hi, top, bed
      real(8) :: d
      character(1024) :: msg

      allocate (plo(nlines), phi(nlines), seen1(nlines), seen2(nlines), flag(nlines))
      allocate (gtop(nlines), gbed(nlines), delta(nlines), zbed(nlines), zc0(nlines))
      plo = 0
      phi = 0
      seen1 = 0
      seen2 = 0

      !...Local search
      do l = 1, nlines
         do i = 1, NVEL
            if (LBCODEI(i) /= 64) cycle
            n1 = NBV(i)
            if (n1 <= 0) cycle
            d = hypot(X(n1) - xq(l), Y(n1) - yq(l))
            if (d > rq(l)) cycle
            g1 = global_node(n1)
            call add_seen(l, g1)
            n2 = IBCONN(i)
            if (n2 <= 0) cycle
            g2 = global_node(n2)
            lo = min(g1, g2)
            hi = max(g1, g2)
            if (plo(l) == 0) then
               plo(l) = lo
               phi(l) = hi
            elseif (plo(l) /= lo .or. phi(l) /= hi) then
               write (msg, '(A,I0,A,4(I0,A))') "TVV: entry ", line_no(l), &
                  ": ambiguous, more than one VEW pair within the search radius (pairs ", &
                  plo(l), "/", phi(l), " and ", lo, "/", hi, ", fulldomain node numbers)."
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
            end if
         end do
      end do

      !...Agreement across ranks: the same pair everywhere it was found
      mlo = plo
      mhi = phi
      call reduce_imax(mlo)
      call reduce_imax(mhi)
      flag = merge(-plo, -huge(0), plo > 0)
      call reduce_imax(flag)
      do l = 1, nlines
         if (mlo(l) == 0) then
            write (msg, '(A,I0,A,ES12.5,A)') "TVV: entry ", line_no(l), &
               ": no IBTYPE=64 VEW pair within the search radius (", rq(l), " m)."
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
         end if
         if (-flag(l) /= mlo(l)) then
            write (msg, '(A,I0,A)') "TVV: entry ", line_no(l), &
               ": ambiguous, subdomains found different VEW pairs within the search radius."
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
         end if
      end do
      flag = merge(-phi, -huge(0), phi > 0)
      call reduce_imax(flag)
      do l = 1, nlines
         if (-flag(l) /= mhi(l)) then
            write (msg, '(A,I0,A)') "TVV: entry ", line_no(l), &
               ": ambiguous, subdomains found different VEW pairs within the search radius."
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
         end if
      end do

      !...Every node seen within the radius must belong to the agreed pair
      flag = 0
      do l = 1, nlines
         if (seen1(l) /= 0 .and. seen1(l) /= mlo(l) .and. seen1(l) /= mhi(l)) flag(l) = 1
         if (seen2(l) /= 0 .and. seen2(l) /= mlo(l) .and. seen2(l) /= mhi(l)) flag(l) = 1
      end do
      call reduce_imax(flag)
      do l = 1, nlines
         if (flag(l) /= 0) then
            write (msg, '(A,I0,A)') "TVV: entry ", line_no(l), &
               ": ambiguous, nodes of another VEW pair are within the search radius."
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
         end if
      end do

      !...Wall-top node (shallower), delta and bed elevation, from a rank holding the pair
      gtop = 0
      delta = -huge(1d0)
      zbed = -huge(1d0)
      zc0 = -huge(1d0)
      do l = 1, nlines
         n1 = local_node(mlo(l))
         n2 = local_node(mhi(l))
         if (n1 <= 0 .or. n2 <= 0) cycle
         i = find_entry(n1, n2)
         if (i == 0) cycle
         if (DP(n1) <= DP(n2)) then
            top = n1
            bed = n2
         else
            top = n2
            bed = n1
         end if
         gtop(l) = global_node(top)
         delta(l) = BARINHT(i) - (-DP(top))
         zbed(l) = -DP(bed)
         zc0(l) = BARINHT(i)
      end do
      call reduce_imax(gtop)
      call reduce_rmax(delta)
      call reduce_rmax(zbed)
      call reduce_rmax(zc0)
      do l = 1, nlines
         if (gtop(l) == 0) then
            write (msg, '(A,I0,A)') "TVV: entry ", line_no(l), &
               ": no subdomain holds both nodes of the VEW pair and its boundary entry."
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
         end if
         gbed(l) = merge(mhi(l), mlo(l), gtop(l) == mlo(l))
      end do

   contains

      subroutine add_seen(l, g)
         integer, intent(in) :: l, g
         if (seen1(l) == 0 .or. seen1(l) == g) then
            seen1(l) = g
         elseif (seen2(l) == 0 .or. seen2(l) == g) then
            seen2(l) = g
         else
            write (msg, '(A,I0,A)') "TVV: entry ", line_no(l), &
               ": ambiguous, nodes of more than one VEW pair are within the search radius."
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
         end if
      end subroutine add_seen

   end subroutine locate_pairs

   !-----------------------------------------------------------------------
   !  Fill the gates with their pairs and the local node and boundary-entry
   !  numbers of each pair (0 where this subdomain does not have them).
   !-----------------------------------------------------------------------
   subroutine build_gates(nlines, line_gate, line_no, gtop, gbed, delta, zbed, zc0)
      integer, intent(in) :: nlines
      integer, intent(in) :: line_gate(:), line_no(:), gtop(:), gbed(:)
      real(8), intent(in) :: delta(:), zbed(:), zc0(:)

      integer :: ig, l, k, m
      character(1024) :: msg

      !...The same pair must not be named twice
      do l = 1, nlines
         do m = l + 1, nlines
            if (gtop(l) == gtop(m)) then
               write (msg, '(A,I0,A,I0,A)') "TVV: entries ", line_no(l), " and ", line_no(m), &
                  " locate the same VEW pair."
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
            end if
         end do
      end do

      do ig = 1, n_tvv_gates
         associate (g => tvv_gates(ig))
            g%npairs = count(line_gate(1:nlines) == ig)
            allocate (g%gtop(g%npairs), g%gbed(g%npairs), g%ltop(g%npairs), g%lbed(g%npairs))
            allocate (g%etop(g%npairs), g%ebed(g%npairs), g%delta(g%npairs), g%zmin(g%npairs))
            allocate (g%zc0(g%npairs))
            k = 0
            do l = 1, nlines
               if (line_gate(l) /= ig) cycle
               k = k + 1
               g%gtop(k) = gtop(l)
               g%gbed(k) = gbed(l)
               g%delta(k) = delta(l)
               g%zmin(k) = zbed(l) + delta(l)
               g%zc0(k) = zc0(l)
               g%ltop(k) = max(local_node(gtop(l)), 0)
               g%lbed(k) = max(local_node(gbed(l)), 0)
               g%etop(k) = 0
               g%ebed(k) = 0
               if (g%ltop(k) > 0 .and. g%lbed(k) > 0) then
                  g%etop(k) = find_entry(g%ltop(k), g%lbed(k))
                  g%ebed(k) = find_entry(g%lbed(k), g%ltop(k))
               end if
            end do
         end associate
      end do

   end subroutine build_gates

   !-----------------------------------------------------------------------
   !  Read the crest table of each gate: a count line, then "time z_c" lines
   !  with strictly increasing times (s). Read from the fulldomain directory.
   !-----------------------------------------------------------------------
   subroutine read_tables()
      use SIZES, only: GLOBALDIR
      use mod_io, only: openFileForRead
      integer :: ig, k, lun, ios
      character(2048) :: fname
      character(1024) :: msg

      lun = 9143
      do ig = 1, n_tvv_gates
         associate (g => tvv_gates(ig))
            fname = trim(GLOBALDIR)//'/'//trim(g%table_file)
            call openFileForRead(lun, trim(fname), ios, required=.true.)
            read (lun, *, iostat=ios) g%nrec
            if (ios /= 0 .or. g%nrec < 1) then
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                              message="TVV: crest table "//trim(fname)// &
                              ": the first line must be the number of records (>= 1).")
            end if
            allocate (g%t(g%nrec), g%z(g%nrec))
            do k = 1, g%nrec
               read (lun, *, iostat=ios) g%t(k), g%z(k)
               if (ios /= 0) then
                  write (msg, '(A,I0)') "TVV: crest table "//trim(fname)//": cannot read record ", k
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
               end if
               if (k > 1) then
                  if (g%t(k) <= g%t(k - 1)) then
                     write (msg, '(A,I0)') "TVV: crest table "//trim(fname)// &
                        ": times must increase strictly; record ", k
                     call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
                  end if
               end if
            end do
            close (lun)
         end associate
      end do

   end subroutine read_tables

   !-----------------------------------------------------------------------
   !  Startup consistency checks (design note, Decision 1-3; input-format
   !  note 4.2). Stops the run on an inconsistency.
   !-----------------------------------------------------------------------
   subroutine check_startup()
      use GLOBAL, only: ILUMP, IHOT, STATIM, TVV_DELTA
      use MESH, only: NP, DP
      use NodalAttributes, only: Tau0, LoadTau0, Tau0DefVal, LoadCondensedNodes, &
                                 NListCondensedNodes, NNodesListCondensedNodes, ListCondensedNodes
      integer :: ig, k, m, n, n1
      integer, allocatable :: topgate(:), toppair(:)
      real(8) :: zc
      character(1024) :: msg

      !...The GWCE left-hand side and TAU0 must not depend on DP
      if (ILUMP == 0) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, message="TVV: time-varying crest VEWs "// &
                        "need a lumped GWCE mass matrix (ILump=1); a consistent mass matrix is not supported yet.")
      end if
      if ((.not. LoadTau0 .and. Tau0 < 0d0) .or. (LoadTau0 .and. Tau0DefVal < 0d0)) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, message="TVV: time-varying crest VEWs "// &
                        "need a constant or prescribed TAU0; depth-dependent or time-varying TAU0 (TAU0 < 0) "// &
                        "is not supported yet.")
      end if

      do ig = 1, n_tvv_gates
         associate (g => tvv_gates(ig))
            do k = 1, g%npairs
               write (msg, '(A,I0,A,I0,A,I0,A)') "TVV: gate ", ig, ", pair ", g%gtop(k), "/", g%gbed(k), ":"
               !...delta: positive, and equal to TVV_DELTA when it is given
               if (is_null(TVV_DELTA)) then
                  if (g%delta(k) <= 0d0) then
                     write (msg, '(A,ES12.5,A)') trim(msg)//" crest minus wall-top elevation is ", &
                        g%delta(k), " m; it must be positive."
                     call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
                  end if
               elseif (abs(g%delta(k) - TVV_DELTA) > TVV_TOL) then
                  write (msg, '(A,ES12.5,A,ES12.5,A)') trim(msg)//" crest minus wall-top elevation is ", &
                     g%delta(k), " m, but TVV_DELTA = ", TVV_DELTA, " m."
                  call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
               end if
               !...Cold start: the crest table must start at the fort.14 crest
               if (IHOT == 0) then
                  zc = tvv_crest_at(ig, STATIM*86400d0)
                  if (abs(zc - g%zc0(k)) > TVV_TOL) then
                     write (msg, '(A,ES12.5,A,ES12.5,A)') trim(msg)//" the crest table gives ", zc, &
                        " m at the start of the run, but the fort.14 crest is ", g%zc0(k), " m."
                     call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
                  end if
               end if
            end do
            !...Table values below a pair's lowest crest (bed + delta) will be clamped
            if (minval(g%z) < maxval(g%zmin) - TVV_TOL) then
               write (msg, '(A,I0,A,ES12.5,A,ES12.5,A)') "TVV: gate ", ig, ": the crest table goes down to ", &
                  minval(g%z), " m, below the lowest crest (bed + delta) of some pairs (up to ", &
                  maxval(g%zmin), " m); the crest will be clamped there."
               call allMessage(WARNING, trim(msg))
            end if
         end associate
      end do

      !...Condensed node groups that contain a wall-top node must consist of
      !   wall-top nodes of one gate with the same DP and delta, so that they
      !   stay level as the gate moves
      if (.not. LoadCondensedNodes) return
      allocate (topgate(NP), toppair(NP))
      topgate = 0
      toppair = 0
      do ig = 1, n_tvv_gates
         do k = 1, tvv_gates(ig)%npairs
            n = tvv_gates(ig)%ltop(k)
            if (n > 0) then
               topgate(n) = ig
               toppair(n) = k
            end if
         end do
      end do
      do k = 1, NListCondensedNodes
         n1 = ListCondensedNodes(k, 1)
         if (all(topgate(ListCondensedNodes(k, 1:NNodesListCondensedNodes(k))) == 0)) cycle
         do m = 1, NNodesListCondensedNodes(k)
            n = ListCondensedNodes(k, m)
            if (topgate(n) == 0 .or. topgate(n) /= topgate(n1) .or. &
                abs(DP(n) - DP(n1)) > TVV_TOL) then
               write (msg, '(A,I0,A)') "TVV: condensed node group ", k, " mixes wall-top nodes of a "// &
                  "time-varying crest gate with other nodes, or its nodes have different depths."
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
            end if
            if (abs(tvv_gates(topgate(n))%delta(toppair(n)) - &
                    tvv_gates(topgate(n1))%delta(toppair(n1))) > TVV_TOL) then
               write (msg, '(A,I0,A)') "TVV: condensed node group ", k, &
                  " has wall-top nodes with different crest offsets (delta)."
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(msg))
            end if
         end do
      end do

   end subroutine check_startup

   !-----------------------------------------------------------------------
   subroutine log_gates()
      use SIZES, only: MYPROC
      integer :: ig, k
      character(4096) :: msg
      character(32) :: buf

      if (MYPROC /= 0) return
      do ig = 1, n_tvv_gates
         associate (g => tvv_gates(ig))
            write (msg, '(A,I0,A,I0,A,2(ES12.5,A),2(ES12.5,A))') "TVV: gate ", ig, ": ", g%npairs, &
               " VEW pairs, crest table '"//trim(g%table_file)//"', crest range [", &
               minval(g%z), ", ", maxval(g%z), "] m, delta range [", &
               minval(g%delta), ", ", maxval(g%delta), "] m"
            call allMessage(INFO, trim(msg))
            msg = "TVV: gate wall-top/bed nodes (fulldomain):"
            do k = 1, g%npairs
               write (buf, '(1X,I0,A,I0)') g%gtop(k), "/", g%gbed(k)
               msg = trim(msg)//trim(buf)
            end do
            call allMessage(INFO, trim(msg))
         end associate
      end do

   end subroutine log_gates

   !-----------------------------------------------------------------------
   !  Local IBTYPE=64 boundary entry with NBV = n1 and IBCONN = n2 (0 if none)
   !-----------------------------------------------------------------------
   integer function find_entry(n1, n2) result(e)
      use BOUNDARIES, only: NVEL, LBCODEI, NBV, IBCONN
      integer, intent(in) :: n1, n2
      integer :: i
      e = 0
      do i = 1, NVEL
         if (LBCODEI(i) == 64 .and. NBV(i) == n1 .and. IBCONN(i) == n2) then
            e = i
            return
         end if
      end do
   end function find_entry

   !-----------------------------------------------------------------------
   !  Fulldomain node number of local node n
   !-----------------------------------------------------------------------
   integer function global_node(n) result(g)
#ifdef CMPI
      use GLOBAL, only: NODES_LG
#endif
      integer, intent(in) :: n
#ifdef CMPI
      g = abs(NODES_LG(n))
#else
      g = n
#endif
   end function global_node

   !-----------------------------------------------------------------------
   !  Local node number of fulldomain node g (0 if not in this subdomain)
   !-----------------------------------------------------------------------
   integer function local_node(g) result(n)
      use MESH, only: NP
#ifdef CMPI
      use GLOBAL, only: NODES_LG
      integer :: i
#endif
      integer, intent(in) :: g
#ifdef CMPI
      n = 0
      do i = 1, NP
         if (abs(NODES_LG(i)) == g) then
            n = i
            return
         end if
      end do
#else
      n = merge(g, 0, g >= 1 .and. g <= NP)
#endif
   end function local_node

   !-----------------------------------------------------------------------
   !  Convert a fort.142 coordinate to model coordinates, as
   !  FIND_BOUNDARY_NODES does (CPP projection when ICS /= 1)
   !-----------------------------------------------------------------------
   subroutine to_model_xy(xin, yin, xout, yout)
      use MESH, only: ICS, SLAM0, SFEA0, DRVSPCOORSROTS, CYLINDERMAP
      use GLOBAL, only: DEG2RAD, IFSPROTS
      real(8), intent(in) :: xin, yin
      real(8), intent(out) :: xout, yout
      real(8) :: lat_temp, lon_temp, latr, lonr

      if (ICS == 1) then
         xout = xin
         yout = yin
      else
         lat_temp = yin*DEG2RAD
         lon_temp = xin*DEG2RAD
         if (IFSPROTS == 1) then
            call DRVSPCOORSROTS(lonr, latr, lon_temp, lat_temp)
         else
            latr = lat_temp
            lonr = lon_temp
         end if
         call CYLINDERMAP(xout, yout, lonr, latr, SLAM0, SFEA0, ICS)
      end if
   end subroutine to_model_xy

   !-----------------------------------------------------------------------
   subroutine nullify_nml()
      X1 = TVV_NULL
      X2 = TVV_NULL
      Y1 = TVV_NULL
      Y2 = TVV_NULL
      ZF = TVV_NULL
      ETA_MAX = TVV_NULL
      TimeStartDay = TVV_NULL
      TimeStartHour = TVV_NULL
      TimeStartMin = TVV_NULL
      TimeStartSec = TVV_NULL
      TimeEndDay = TVV_NULL
      TimeEndHour = TVV_NULL
      TimeEndMin = TVV_NULL
      TimeEndSec = TVV_NULL
      FailureDurationDay = TVV_NULL
      FailureDurationHour = TVV_NULL
      FailureDurationMin = TVV_NULL
      FailureDurationSec = TVV_NULL
      SearchRadius = TVV_NULL
      VaryType = -99999
      HOT = -99999
      LOOP = -99999
      NLOOPS = -99999
      ScheduleFile = "NOFILE"
   end subroutine nullify_nml

   pure logical function is_null(r)
      real(8), intent(in) :: r
      is_null = abs(r - TVV_NULL) <= epsilon(1d0)
   end function is_null

   !-----------------------------------------------------------------------
   !  In-place maximum over all ranks (no-op in serial)
   !-----------------------------------------------------------------------
   subroutine reduce_imax(a)
#ifdef CMPI
      use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_INTEGER, MPI_MAX
      use GLOBAL, only: COMM
#endif
      integer, intent(inout) :: a(:)
#ifdef CMPI
      if (size(a) > 0) call MPI_Allreduce(MPI_IN_PLACE, a, size(a), MPI_INTEGER, MPI_MAX, COMM)
#endif
   end subroutine reduce_imax

   subroutine reduce_rmax(a)
#ifdef CMPI
      use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX
      use GLOBAL, only: COMM
#endif
      real(8), intent(inout) :: a(:)
#ifdef CMPI
      if (size(a) > 0) call MPI_Allreduce(MPI_IN_PLACE, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, COMM)
#endif
   end subroutine reduce_rmax

end module mod_tvv
