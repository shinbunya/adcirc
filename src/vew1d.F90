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
!  MODULE VEW1D
!
!  This module declares subroutines for vertical element walls and
!  1D condensation
!
!  2022-09-21: Created by SB
!
!*******************************************************************
module vew1d
   implicit none

   private

   public :: ROTATE_AT_CONDENSEDNODES_ALL, ROTATEBACK_AT_CONDENSEDNODES_ALL, &
             REMOVE_NORMAL_AT_CONDENSEDNODES, GET_STREAMDIRECTION_AT_ELEMENT, &
             COMPUTE_VEW_ADVECT_VELOCITY, UADV, VADV, &
             VEW_ADVECT_FR_ATTACH, VEW_ADVECT_FR_FULL, VEW1D_ROTATE_VELOCITY, &
             SETUP_COUPLED_GROUPS, UPDATE_ACTIVE_GROUPS, GROUP_LOAD_SUM, GROUP_COPY_INCREMENT, &
             GROUP_SLOT, NSGROUPS, NSGROUPS_GBL, ISSGMEMBER, SGELEM, GMAP, NeiTabG, NNeighG, &
             GROUP_BLOCK_RESET, GROUP_BLOCK_ADD, GROUP_IMAX

   !...VEW1D_ROTATE_VELOCITY (VEW1DChannelControl, default off): rotate the element velocities and
   !   fluxes at two-node condensed groups along uncoupled VEW or land boundaries to the stream
   !   direction in the GWCE and momentum element integrals (ROTATE_AT_CONDENSEDNODES). The nodal
   !   treatments (rotate-back and removal of the normal component of the momentum load, wind stress
   !   and Coriolis) do not depend on it.
   logical :: VEW1D_ROTATE_VELOCITY = .false.

   !...Coupled groups for the consistent-mass GWCE (ILump = 0). Nodes that share one water level
   !   (VEW pairs in Case 1 with both nodes wet, wet members of a condensed node group, and chains of
   !   these) get one unknown: every element matrix is assembled at the rows and columns of the
   !   group representatives (GMAP), which gives P^T K P, the members' loads are summed to the
   !   representative and its increment is copied back to the members after JCG.
   !   Static groups: connected sets under all IBTYPE=64 pairs and all condensed groups, i.e. every
   !   coupling that can become active. NeiTabG/NNeighG is the neighbour table of the GWCE matrix:
   !   NeiTab plus, for each static member, the other members and all members' neighbours, and for
   !   each of those neighbours, all members. Without static groups it is a copy of NeiTab.
   !   Formulation: docs/reports/vew-explicit-and-implicit-formulation (operable-structure project).
   integer :: NSGROUPS = 0                      ! number of static groups
   integer :: NSGROUPS_GBL = 0                  ! its maximum over the ranks (set by initGWCE)
   logical, allocatable :: ISSGMEMBER(:)        ! node belongs to a static group
   logical, allocatable :: SGELEM(:)            ! element has a node in a static group
   integer, allocatable :: GMAP(:)              ! representative of the node's active group
   integer, allocatable :: GMAPPREV(:)          ! GMAP of the previous step
   integer, allocatable :: SGMEMBERS(:)         ! static members
   integer, allocatable :: NeiTabG(:, :), NNeighG(:)
   integer, allocatable :: GPARENT(:)           ! union-find work array
   !...entries of the GWCE matrix changed by GROUP_BLOCK_ADD since the last reset (row, column
   !   slot, value before the change), so that the assembled matrix can be restored
   integer :: NBLK = 0
   integer, allocatable :: BLKROW(:), BLKSLOT(:)
   real(8), allocatable :: BLKORIG(:)
   real(8), allocatable :: GSUM(:)              ! work array for the group sums
   integer, allocatable :: GCNT(:)

   !...Recovered crest velocity at coupled (submerged) IBTYPE=64 VEW node pairs, on with
   !   VEW_CREST_HYDRAULICS (WEIR module, VEW1DChannelControl namelist). A coupled pair shares one
   !   velocity, so the shallower (wall-top) side carries less than the discharge reaching the wall.
   !   The momentum and GWCE advective terms then use, at the nodes of each coupled pair, a velocity
   !   whose wall-normal component (boundary normal CSII/SIII of the pair) is
   !      v'n = q / min(H1, H2),
   !   q the wall-normal discharge averaged over the wet neighbours of the pair's foot (deeper) node
   !   that are not on an IBTYPE=64 boundary; the tangential component is kept as is. Where the foot
   !   has no such neighbour, or both nodes have the same depth, v'n = vn*0.5*(H1 + H2)/min(H1, H2).
   !   UADV/VADV hold the velocities those terms use. Further:
   !   - Lee loss: the foot node of a pair whose bed side is downstream (flow leaving the wall, judged
   !     from the neighbours' flow) keeps its own velocity. The approach side and the whole crest carry
   !     v', the lee-side bed elements see the v' -> v drop, so the velocity head over the crest is lost
   !     as in the broad-crested weir formulas (expansion loss) instead of recovered.
   !   - Froude taper: with Fr = |v'n|/sqrt(g*min(H)), the recovery has weight 0 (flow attached to the
   !     wall, no recovery and no lee loss) below VEW_ADVECT_FR_ATTACH, 1 above VEW_ADVECT_FR_FULL,
   !     linear between.
   !   - Froude cap: |v'n| <= sqrt(g*min(H1, H2)), so the recovered flow stays subcritical.
   !   - Condensed wall tops: both faces of a wall share one discharge, the one of the face
   !     whose foot is on the upstream side, stored per condensed group (QGROUP, keyed by the group's
   !     first node). While a pair of a wall-top group is in weir mode (uncoupled with both nodes
   !     wet), the coupled pairs of that group get no recovered velocity (GROUPWEIR), so the approach
   !     to the wall is not drawn down onto the crest a second time: the weir formula on the lee pair
   !     already represents the crest.
   real(8), allocatable, target :: UADV(:), VADV(:)
   real(8) :: VEW_ADVECT_FR_ATTACH = 0.05d0
   real(8) :: VEW_ADVECT_FR_FULL = 0.1d0
   real(8), allocatable :: QGROUP(:)
   logical, allocatable :: QGROUPSET(:)
   logical, allocatable :: ISVEWNODE(:)
   logical, allocatable :: GROUPWEIR(:)

contains

   !-------------------------------------------------------------
   subroutine GET_STREAMDIRECTION_AT_ELEMENT &
      (NM1, NM2, NM3, IA1, IA2, IB1, IB2, &
       CAX, CAY, CBX, CBY, SX, SY, STATUS)
      use MESH, only: X, Y, findIBTYPEAtNode
      use NodalAttributes, only: NCondensedNodes, CondensedNodes

      implicit none

      integer, intent(IN) :: NM1, NM2, NM3
      integer, intent(OUT) :: IA1, IA2, IB1, IB2
      real(8), intent(OUT) :: CAX, CAY, CBX, CBY, SX, SY
      integer, intent(OUT) :: STATUS
      integer :: IBTYPE_A1, IBTYPE_A2, IBTYPE_B1
      real(8) :: LEN
      logical :: ON_BOUNDARY

      STATUS = 0

      if ((NCondensedNodes(nm1) < 2 .or. & ! if any of the nodes are not condensed nodes
           NCondensedNodes(nm2) < 2 .or. &
           NCondensedNodes(nm3) < 2) .or. & ! or
          (NCondensedNodes(nm1) > 2 .and. & ! if this is a joint element
           NCondensedNodes(nm2) > 2 .and. &
           NCondensedNodes(nm3) > 2)) then
         status = -1
         return
      end if

      if (CondensedNodePairFound(NM1, NM2)) then
         IA1 = NM1
         IA2 = NM2
         IB1 = NM3
      else if (CondensedNodePairFound(NM2, NM1)) then
         IA1 = NM2
         IA2 = NM1
         IB1 = NM3
      else if (CondensedNodePairFound(NM2, NM3)) then
         IA1 = NM2
         IA2 = NM3
         IB1 = NM1
      else if (CondensedNodePairFound(NM3, NM2)) then
         IA1 = NM3
         IA2 = NM2
         IB1 = NM1
      else if (CondensedNodePairFound(NM3, NM1)) then
         IA1 = NM3
         IA2 = NM1
         IB1 = NM2
      else if (CondensedNodePairFound(NM1, NM3)) then
         IA1 = NM1
         IA2 = NM3
         IB1 = NM2
      else
         status = -1
         return !when condensed node pair not found
      end if
      IB2 = FIND_IB2(IB1, IA1, IA2)
      if (IB2 < 1) then
         status = -1
         return !when condensed node pair not found
      end if
      if (CondensedNodes(IB2, 1) == 0) then
         status = -1
         return !when condensed node pair not found
      end if

      ibtype_a1 = findIBTYPEAtNode(ia1)
      ibtype_a2 = findIBTYPEAtNode(ia2)
      ibtype_b1 = findIBTYPEAtNode(ib1)

      on_boundary = .false.
      if ((ibtype_a1 == 20 .and. ibtype_b1 == 20) .or. &
          (ibtype_a2 == 20 .and. ibtype_b1 == 20) .or. &
          (ibtype_a1 == 64 .and. ibtype_b1 == 64) .or. &
          (ibtype_a2 == 64 .and. ibtype_b1 == 64)) then
         on_boundary = .true.
      end if

      if (.not. on_boundary) then
         status = -1
         return
      end if

      CAX = 0.5d0*(X(IA1) + X(IA2))
      CAY = 0.5d0*(Y(IA1) + Y(IA2))
      CBX = 0.5d0*(X(IB1) + X(IB2))
      CBY = 0.5d0*(Y(IB1) + Y(IB2))
      SX = CAX - CBX
      SY = CAY - CBY
      LEN = sqrt(SX*SX + SY*SY)
      SX = SX/LEN
      SY = SY/LEN

      return
   end subroutine GET_STREAMDIRECTION_AT_ELEMENT

   !-------------------------------------------------------------
   subroutine ROTATE_AT_CONDENSEDNODES(VECX, VECY, IND, SX, SY)
      use NodalAttributes, only: LoadCondensedNodes, NCondensedNodes
      use Mesh, only: LBArray_Pointer
      implicit none

      real(8), intent(INOUT) :: VECX, VECY
      real(8), intent(IN) :: SX, SY
      integer, intent(IN) :: IND
      real(8) :: CS, SI, VLEN, P
      integer :: IBND
      logical :: DO_ROTATE

      if (.not. VEW1D_ROTATE_VELOCITY) return
      if (.not. LoadCondensedNodes) return
      if (NCondensedNodes(IND) /= 2) return

      IBND = LBArray_Pointer(IND)
      if (IBND <= 0) return

      ! Rotate the vector if there is no bank top
      ! side, i.e., along land boundaries, or along not-submerged
      ! vertical element walls
      DO_ROTATE = ROTATE_OR_NOT(IBND)

      if (.not. DO_ROTATE) return

      CS = SX
      SI = SY

      VLEN = sqrt(VECX*VECX + VECY*VECY)

      P = VECX*CS + VECY*SI

      VECX = VLEN*CS*sign(1.d0, P)
      VECY = VLEN*SI*sign(1.d0, P)

   end subroutine ROTATE_AT_CONDENSEDNODES

   !-------------------------------------------------------------
   subroutine ROTATEBACK_AT_CONDENSEDNODES(VECX, VECY, IND)
      use NodalAttributes, only: LoadCondensedNodes, NCondensedNodes
      use Mesh, only: LBArray_Pointer
      use BOUNDARIES, only: CSIICN, SIIICN
      implicit none

      real(8), intent(INOUT) :: VECX, VECY
      integer, intent(IN) :: IND
      real(8) :: CS, SI, VLEN, P
      integer :: IBND
      logical :: DO_ROTATE

      if (.not. LoadCondensedNodes) return
      if (NCondensedNodes(IND) /= 2) return

      IBND = LBArray_Pointer(IND)
      if (IBND <= 0) return

      ! Rotate back the vector if there is no bank top
      ! side, i.e., along land boundaries, or along not-submerged
      ! vertical element walls
      DO_ROTATE = ROTATE_OR_NOT(IBND)

      if (.not. DO_ROTATE) return

      VLEN = sqrt(VECX*VECX + VECY*VECY)
      CS = -SIIICN(IBND)
      SI = CSIICN(IBND)
      P = VECX*CS + VECY*SI
      VECX = VLEN*CS*sign(1.d0, P)
      VECY = VLEN*SI*sign(1.d0, P)

      return
   end subroutine ROTATEBACK_AT_CONDENSEDNODES

   !-------------------------------------------------------------
   subroutine ROTATE_AT_CONDENSEDNODES_ALL &
      (NM1, NM2, NM3, &
       U1N1, V1N1, U1N2, V1N2, U1N3, V1N3, &
       QX1N1, QY1N1, QX1N2, QY1N2, QX1N3, QY1N3, &
       STATUS)
      !-------------------------------------------------------------
      implicit none

      integer, intent(IN) :: NM1, NM2, NM3
      real(8), intent(INOUT) :: &
         U1N1, V1N1, U1N2, V1N2, U1N3, V1N3, &
         QX1N1, QY1N1, QX1N2, QY1N2, QX1N3, QY1N3
      integer, intent(OUT) :: STATUS
      real(8) :: CAX, CAY, CBX, CBY, SX, SY
      integer :: IA1, IA2, IB1, IB2

      call GET_STREAMDIRECTION_AT_ELEMENT &
         (NM1, NM2, NM3, IA1, IA2, IB1, IB2, &
          CAX, CAY, CBX, CBY, SX, SY, STATUS)

      if (STATUS /= 0) return

      call ROTATE_AT_CONDENSEDNODES(U1N1, V1N1, NM1, SX, SY)
      call ROTATE_AT_CONDENSEDNODES(QX1N1, QY1N1, NM1, SX, SY)

      call ROTATE_AT_CONDENSEDNODES(U1N2, V1N2, NM2, SX, SY)
      call ROTATE_AT_CONDENSEDNODES(QX1N2, QY1N2, NM2, SX, SY)

      call ROTATE_AT_CONDENSEDNODES(U1N3, V1N3, NM3, SX, SY)
      call ROTATE_AT_CONDENSEDNODES(QX1N3, QY1N3, NM3, SX, SY)

      return
   end subroutine ROTATE_AT_CONDENSEDNODES_ALL

   !-------------------------------------------------------------
   subroutine ROTATEBACK_AT_CONDENSEDNODES_ALL &
      (NM1, NM2, NM3, TEMP_LV_A1, TEMP_LV_B1, &
       TEMP_LV_A2, TEMP_LV_B2, TEMP_LV_A3, TEMP_LV_B3, &
       STATUS)
      !-------------------------------------------------------------
      implicit none

      integer, intent(IN) :: NM1, NM2, NM3
      real(8), intent(INOUT) :: &
         TEMP_LV_A1, TEMP_LV_B1, TEMP_LV_A2, TEMP_LV_B2, &
         TEMP_LV_A3, TEMP_LV_B3
      integer, intent(OUT) :: STATUS
      real(8) :: CAX, CAY, CBX, CBY, SX, SY
      integer :: IA1, IA2, IB1, IB2

      call GET_STREAMDIRECTION_AT_ELEMENT &
         (NM1, NM2, NM3, IA1, IA2, IB1, IB2, &
          CAX, CAY, CBX, CBY, SX, SY, STATUS)

      if (STATUS /= 0) return

      call ROTATEBACK_AT_CONDENSEDNODES(TEMP_LV_A1, TEMP_LV_B1, NM1)
      call ROTATEBACK_AT_CONDENSEDNODES(TEMP_LV_A2, TEMP_LV_B2, NM2)
      call ROTATEBACK_AT_CONDENSEDNODES(TEMP_LV_A3, TEMP_LV_B3, NM3)

      return
   end subroutine ROTATEBACK_AT_CONDENSEDNODES_ALL

   !-------------------------------------------------------------
   subroutine REMOVE_NORMAL_AT_CONDENSEDNODES(VECX, VECY, IND, STATUS)
      use NodalAttributes, only: LoadCondensedNodes, NCondensedNodes
      use Mesh, only: LBArray_Pointer
      use BOUNDARIES, only: CSIICN, SIIICN

      implicit none

      real(8), intent(INOUT) :: VECX, VECY
      integer, intent(IN) :: IND
      integer, intent(OUT) :: STATUS
      real(8) :: CS, SI, CSCS, CSSI, SISI, fBUF1, fBUF2
      integer :: IBND
      logical :: DO_ROTATE

      STATUS = -1

      if (.not. LoadCondensedNodes) return
      if (NCondensedNodes(IND) /= 2) return

      IBND = LBArray_Pointer(IND)
      if (IBND <= 0) return

      ! Remove the normal component if there is no water on the
      ! bank, i.e., along land boundaries, or along not-submerged
      ! vertical element walls
      DO_ROTATE = ROTATE_OR_NOT(IBND)

      if (.not. DO_ROTATE) return

      CS = -SIIICN(IBND)
      SI = CSIICN(IBND)
      CSCS = CS*CS
      CSSI = CS*SI
      SISI = SI*SI

      fBUF1 = VECX*CSCS + VECY*CSSI
      fBUF2 = VECX*CSSI + VECY*SISI

      VECX = fBUF1
      VECY = fBUF2

      STATUS = 0

      return
   end subroutine REMOVE_NORMAL_AT_CONDENSEDNODES

   !-------------------------------------------------------------
   ! FUNCTION ROTATE_OR_NOT
   !
   ! Judges if velocity solutions should be rotated at a node
   ! along a boundary. They should be rotated if 1) a node is
   ! along a IBTYPE=64 vertical element wall boundary and the
   ! boundary is not submerged, or if 2) a node is along any
   ! zero normal flow boundary. This function is applied to a
   ! node where the condensed_nodes nodal attribute is given.
   !-------------------------------------------------------------
   function ROTATE_OR_NOT(IBND)
      use GLOBAL, only: NODECODE
      use BOUNDARIES, only: IBCONN, LBCODEI, NBV, ISSUBMERGED64

      implicit none

      integer, intent(IN) :: IBND
      integer :: NNBB1, NNBB2
      logical :: ROTATE_OR_NOT

      ROTATE_OR_NOT = .false.
      if (LBCODEI(IBND) == 64) then
         NNBB1 = NBV(IBND)
         NNBB2 = IBCONN(IBND)
         if ((ISSUBMERGED64(IBND) == 0) .or. &
             (NODECODE(NNBB1) == 0 .or. NODECODE(NNBB2) == 0)) then
            ROTATE_OR_NOT = .true.
         end if
      else if (LBCODEI(IBND) == 0 .or. LBCODEI(IBND) == 1 .or. &
               LBCODEI(IBND) == 10 .or. LBCODEI(IBND) == 1 .or. &
               LBCODEI(IBND) == 20 .or. LBCODEI(IBND) == 21) then
         ROTATE_OR_NOT = .true.
      end if
   end function ROTATE_OR_NOT

   function CondensedNodePairFound(N1, N2)
      use NodalAttributes, only: CondensedNodes, CondensedNodesNoOfVals
      implicit none
      integer, intent(in) :: N1, N2
      logical :: CondensedNodePairFound
      integer :: I, N

      CondensedNodePairFound = .false.
      if (CondensedNodes(N1, 1) < 0) then
         N = abs(CondensedNodes(N1, 1))
      else
         N = N1
      end if
      do I = 1, CondensedNodesNoOfVals
         if (CondensedNodes(N, I) == N2) then
            CondensedNodePairFound = .true.
            exit
         end if
      end do

   end function CondensedNodePairFound

   function FIND_IB2(N1, N2, N3)
      use MESH, only: NNEIGHELE, NEITABELE, NM
      implicit none
      integer, intent(IN) :: N1, N2, N3
      integer :: FIND_IB2
      integer :: nnm1, nnm2, nnm3, tmp
      integer :: I, IE
      logical :: FOUND

      find_ib2 = -1
      found = .false.

      do i = 1, nneighele(n1)
         ie = neitabele(n1, i)
         if (ie <= 0) cycle ! sb251003: Added this check to avoid out-of-bounds access.
         nnm1 = nm(ie, 1)
         nnm2 = nm(ie, 2)
         nnm3 = nm(ie, 3)
         if (nnm1 == n1) then
            continue ! do nothing
         else if (nnm2 == n1) then
            tmp = nnm1
            nnm1 = nnm2
            nnm2 = nnm3
            nnm3 = tmp
         else if (nnm3 == n1) then
            tmp = nnm1
            nnm1 = nnm3
            nnm3 = nnm2
            nnm2 = tmp
         else
            stop 'Error in FIND_IB2: node not found in neighbor element'
         end if
         if ((nnm2 == n2 .and. nnm3 == n3) .or. &
             (nnm2 == n3 .and. nnm3 == n2)) cycle
         if (nnm2 == n2 .or. nnm2 == n3) then
            find_ib2 = nnm3
            found = .true.
         else if (nnm3 == n2 .or. nnm3 == n3) then
            find_ib2 = nnm2
            found = .true.
         end if
         if (found) exit
      end do
   end function FIND_IB2

   !-------------------------------------------------------------
   !  COMPUTE_VEW_ADVECT_VELOCITY
   !  Fill UADV/VADV with UU1/VV1 and, if VEW_CREST_HYDRAULICS is set, overwrite them
   !  at both nodes of each coupled IBTYPE=64 pair (both nodes wet, min depth > H0)
   !  with the recovered velocity above. Called once per time step before the GWCE.
   !-------------------------------------------------------------
   subroutine COMPUTE_VEW_ADVECT_VELOCITY()
      use GLOBAL, only: UU1, VV1, ETA2, NODECODE, IFNLFA, H0
      use ADC_CONSTANTS, only: G
      use MESH, only: NP, DP, X, Y, NNeigh, NeiTab
      use BOUNDARIES, only: NBOU, NVELL, LBCODEI, NBV, IBCONN, ISSUBMERGED64, CSII, SIII
      use NodalAttributes, only: LoadCondensedNodes, CondensedNodes
      use WEIR, only: VEW_CREST_HYDRAULICS
      implicit none
      integer :: I, J, K, N1, N2, NF, JJ, NT, KG
      real(8) :: HA, HB, HMIN, FAC, UAVG, VAVG, UN, U0, V0, DX, DY, QN, VNP, VCRIT
      real(8) :: UNB, VNB, WREC, FRC
      integer :: NB, CNT, II
      logical :: LEEFOOT

      if (.not. allocated(UADV)) allocate (UADV(NP), VADV(NP))
      UADV(:) = UU1(:)
      VADV(:) = VV1(:)
      if (.not. VEW_CREST_HYDRAULICS) return

      if (.not. allocated(ISVEWNODE)) then
         allocate (ISVEWNODE(NP))
         ISVEWNODE(:) = .false.
         II = 0
         do K = 1, NBOU
            select case (LBCODEI(II + 1))
            case (64)
               do J = 1, NVELL(K)
                  II = II + 1
                  ISVEWNODE(NBV(II)) = .true.
                  ISVEWNODE(IBCONN(II)) = .true.
               end do
               II = II + NVELL(K)
            case (4, 24, 5, 25)
               II = II + NVELL(K)*2
            case default
               II = II + NVELL(K)
            end select
         end do
      end if

      !...Pass 1 (condensed wall tops): discharge of the upstream face of each condensed wall-top group, and the
      !   groups with a pair in weir mode
      if (LoadCondensedNodes) then
         if (.not. allocated(QGROUP)) allocate (QGROUP(NP), QGROUPSET(NP), GROUPWEIR(NP))
         QGROUPSET(:) = .false.
         GROUPWEIR(:) = .false.
         I = 0
         do K = 1, NBOU
            select case (LBCODEI(I + 1))
            case (64)
               do J = 1, NVELL(K)
                  I = I + 1
                  N1 = NBV(I)
                  N2 = IBCONN(I)
                  if (DP(N1) == DP(N2)) cycle
                  if (NODECODE(N1) == 0 .or. NODECODE(N2) == 0) cycle
                  if (ISSUBMERGED64(I) == 0) then
                     KG = CONDENSED_GROUP(merge(N1, N2, DP(N2) > DP(N1)))
                     if (KG > 0) GROUPWEIR(KG) = .true.
                     cycle
                  end if
                  NF = merge(N2, N1, DP(N2) > DP(N1))
                  NT = merge(N1, N2, NF == N2)
                  KG = CONDENSED_GROUP(NT)
                  if (KG == 0) cycle
                  call NEIGHBOR_FLOW(NF, I, QN, UNB, VNB, DX, DY, CNT)
                  if (CNT == 0) cycle
                  !...upstream face: the flow goes from the foot's neighbours towards the wall
                  if (UNB*DX + VNB*DY < 0.d0) then
                     QGROUP(KG) = QN
                     QGROUPSET(KG) = .true.
                  end if
               end do
               I = I + NVELL(K)
            case (4, 24, 5, 25)
               I = I + NVELL(K)*2
            case default
               I = I + NVELL(K)
            end select
         end do
      end if

      I = 0
      do K = 1, NBOU
         select case (LBCODEI(I + 1))
         case (64)
            do J = 1, NVELL(K)
               I = I + 1
               if (ISSUBMERGED64(I) == 0) cycle
               N1 = NBV(I)
               N2 = IBCONN(I)
               if (NODECODE(N1) == 0 .or. NODECODE(N2) == 0) cycle
               HA = DP(N1) + IFNLFA*ETA2(N1)
               HB = DP(N2) + IFNLFA*ETA2(N2)
               HMIN = min(HA, HB)
               if (HMIN <= H0) cycle
               if (LoadCondensedNodes .and. DP(N1) /= DP(N2)) then
                  KG = CONDENSED_GROUP(merge(N1, N2, DP(N2) > DP(N1)))
                  if (KG > 0) then
                     if (GROUPWEIR(KG)) cycle ! plain velocities at both nodes
                  end if
               end if
               FAC = 0.5d0*(HA + HB)/HMIN
               UAVG = 0.5d0*(UU1(N1) + UU1(N2))
               VAVG = 0.5d0*(VV1(N1) + VV1(N2))
               U0 = UAVG
               V0 = VAVG
               !...Lee-side foot: the deeper node of the pair, if the flow goes from it into the
               !   bed-side domain. Its neighbours all lie on the bed side (the coincident wall-top
               !   node is not connected to it), so their mean position minus its own points into it.
               LEEFOOT = .false.
               if (DP(N1) /= DP(N2)) then
                  NF = N1
                  if (DP(N2) > DP(N1)) NF = N2
                  DX = 0.d0
                  DY = 0.d0
                  do JJ = 2, NNeigh(NF)
                     DX = DX + X(NeiTab(NF, JJ)) - X(NF)
                     DY = DY + Y(NeiTab(NF, JJ)) - Y(NF)
                  end do
                  LEEFOOT = (U0*DX + V0*DY) > 0.d0
               end if
               !...scale only the normal component: u' = u + (FAC - 1)*(u.n)n
               UN = UAVG*CSII(I) + VAVG*SIII(I)
               UAVG = UAVG + (FAC - 1.d0)*UN*CSII(I)
               VAVG = VAVG + (FAC - 1.d0)*UN*SIII(I)
               if (DP(N1) /= DP(N2)) then
                  NF = N1
                  if (DP(N2) > DP(N1)) NF = N2
                  call NEIGHBOR_FLOW(NF, I, QN, UNB, VNB, DX, DY, CNT)
                  if (CNT > 0) then
                     !...one discharge per wall with condensed tops: the upstream face's (magnitude), with
                     !   the sign of this face's normal flux
                     if (LoadCondensedNodes) then
                        KG = CONDENSED_GROUP(merge(N1, N2, NF == N2))
                        if (KG > 0) then
                           if (QGROUPSET(KG)) QN = sign(abs(QGROUP(KG)), QN)
                        end if
                     end if
                     VNP = QN/HMIN
                     UN = U0*CSII(I) + V0*SIII(I)
                     UAVG = U0 + (VNP - UN)*CSII(I)
                     VAVG = V0 + (VNP - UN)*SIII(I)
                     !...lee side from the neighbours' flow, not the pair's own velocity
                     LEEFOOT = (UNB*DX + VNB*DY) > 0.d0
                     !...taper with the crest Froude number
                     FRC = abs(VNP)/sqrt(G*HMIN)
                     if (FRC <= VEW_ADVECT_FR_ATTACH) then
                        WREC = 0.d0
                     elseif (FRC >= VEW_ADVECT_FR_FULL) then
                        WREC = 1.d0
                     else
                        WREC = (FRC - VEW_ADVECT_FR_ATTACH)/(VEW_ADVECT_FR_FULL - VEW_ADVECT_FR_ATTACH)
                     end if
                     if (WREC <= 0.d0) cycle ! attached: plain velocities at both nodes, no lee loss
                     UAVG = U0 + WREC*(UAVG - U0)
                     VAVG = V0 + WREC*(VAVG - V0)
                  end if
               end if
               !...Froude cap
               VCRIT = sqrt(G*HMIN)
               UN = UAVG*CSII(I) + VAVG*SIII(I)
               if (abs(UN) > VCRIT) then
                  VNP = sign(VCRIT, UN)
                  UAVG = UAVG + (VNP - UN)*CSII(I)
                  VAVG = VAVG + (VNP - UN)*SIII(I)
               end if
               if (LEEFOOT) then
                  !...recovered velocity at the wall-top node only; the lee foot keeps its own
                  N1 = merge(N2, N1, NF == N1) ! N1 := the wall-top node
                  UADV(N1) = UAVG
                  VADV(N1) = VAVG
               else
                  UADV(N1) = UAVG
                  UADV(N2) = UAVG
                  VADV(N1) = VAVG
                  VADV(N2) = VAVG
               end if
            end do
            I = I + NVELL(K)
         case (4, 24, 5, 25)
            I = I + NVELL(K)*2
         case default
            I = I + NVELL(K)
         end select
      end do
   end subroutine COMPUTE_VEW_ADVECT_VELOCITY

   !-------------------------------------------------------------
   !  NEIGHBOR_FLOW: mean wall-normal discharge QN (boundary normal of entry IB) and mean velocity
   !  (UNB, VNB) over the wet neighbours of the foot node NF that are not on a VEW (CNT of them),
   !  and the sum (DX, DY) of the offsets of all its neighbours (points from NF into its bed side).
   !-------------------------------------------------------------
   subroutine NEIGHBOR_FLOW(NF, IB, QN, UNB, VNB, DX, DY, CNT)
      use GLOBAL, only: UU1, VV1, ETA2, NODECODE, IFNLFA
      use MESH, only: DP, X, Y, NNeigh, NeiTab
      use BOUNDARIES, only: CSII, SIII
      implicit none
      integer, intent(IN) :: NF, IB
      real(8), intent(OUT) :: QN, UNB, VNB, DX, DY
      integer, intent(OUT) :: CNT
      integer :: JJ, NB

      QN = 0.d0
      UNB = 0.d0
      VNB = 0.d0
      DX = 0.d0
      DY = 0.d0
      CNT = 0
      do JJ = 2, NNeigh(NF)
         NB = NeiTab(NF, JJ)
         DX = DX + X(NB) - X(NF)
         DY = DY + Y(NB) - Y(NF)
         if (NB == NF .or. ISVEWNODE(NB) .or. NODECODE(NB) == 0) cycle
         QN = QN + (UU1(NB)*CSII(IB) + VV1(NB)*SIII(IB))*(DP(NB) + IFNLFA*ETA2(NB))
         UNB = UNB + UU1(NB)
         VNB = VNB + VV1(NB)
         CNT = CNT + 1
      end do
      if (CNT > 0) then
         QN = QN/CNT
         UNB = UNB/CNT
         VNB = VNB/CNT
      end if
   end subroutine NEIGHBOR_FLOW

   !-------------------------------------------------------------
   !  CONDENSED_GROUP: first node of the condensed node group of node N (0 if not condensed)
   !-------------------------------------------------------------

   integer function CONDENSED_GROUP(N)
      use NodalAttributes, only: CondensedNodes
      implicit none
      integer, intent(IN) :: N
      CONDENSED_GROUP = 0
      if (CondensedNodes(N, 1) > 0) then
         CONDENSED_GROUP = N
      elseif (CondensedNodes(N, 1) < 0) then
         CONDENSED_GROUP = -CondensedNodes(N, 1)
      end if
   end function CONDENSED_GROUP

   !-------------------------------------------------------------
   !  SETUP_COUPLED_GROUPS: static groups and the GWCE neighbour table for the consistent-mass GWCE.
   !  Called once from READ_INPUT (ILump = 0), after the mesh, the boundaries and the condensed node
   !  groups are set up and before the GWCE matrix is allocated. Sets MNEIG, the number of columns of
   !  COEF.
   !-------------------------------------------------------------
   subroutine SETUP_COUPLED_GROUPS()
      use SIZES, only: MNP, MNEI, MNEIG
      use MESH, only: NP, NE, NM, NeiTab, NNeigh
      use BOUNDARIES, only: NVEL, NBV, IBCONN, LBCODEI
      use NodalAttributes, only: LoadCondensedNodes, NListCondensedNodes, &
                                 NNodesListCondensedNodes, ListCondensedNodes
      use WEIR, only: WEIR_SEMI_IMPLICIT
      use mod_logging, only: allMessage, INFO
      implicit none
      integer :: I, J, K, L, N, R, NMEM, NCOL, IE, NWP
      integer, allocatable :: ROOT(:), SIZEOF(:), FIRST(:), NEXT(:)
      integer, allocatable :: TAB(:, :), CNT(:), UNION(:)
      integer :: NU, NU0
      character(len=256) :: MSG

      allocate (GPARENT(NP), GMAP(NP), GMAPPREV(NP), ISSGMEMBER(NP), GSUM(NP), GCNT(NP))
      allocate (SGELEM(NE))
      do I = 1, NP
         GPARENT(I) = I
         GMAP(I) = I
      end do
      GMAPPREV = 0
      ISSGMEMBER = .false.
      SGELEM = .false.

      !...static groups: all IBTYPE=64 pairs and all condensed groups
      do I = 1, NVEL
         if (LBCODEI(I) /= 64) cycle
         if (NBV(I) <= 0 .or. IBCONN(I) <= 0) cycle
         call GJOIN(NBV(I), IBCONN(I))
      end do
      if (LoadCondensedNodes) then
         do K = 1, NListCondensedNodes
            if (ListCondensedNodes(K, 1) == 0) cycle
            do L = 2, NNodesListCondensedNodes(K)
               call GJOIN(ListCondensedNodes(K, 1), ListCondensedNodes(K, L))
            end do
         end do
      end if

      allocate (ROOT(NP), SIZEOF(NP), FIRST(NP), NEXT(NP))
      SIZEOF = 0
      do I = 1, NP
         ROOT(I) = GFIND(I)
         SIZEOF(ROOT(I)) = SIZEOF(ROOT(I)) + 1
      end do
      !...member lists per group, in increasing node order (FIRST(root) -> NEXT)
      FIRST = 0
      NEXT = 0
      NMEM = 0
      NSGROUPS = 0
      do I = NP, 1, -1
         R = ROOT(I)
         if (SIZEOF(R) < 2) cycle
         NEXT(I) = FIRST(R)
         FIRST(R) = I
         ISSGMEMBER(I) = .true.
         NMEM = NMEM + 1
      end do
      do I = 1, NP
         if (SIZEOF(I) >= 2 .and. ROOT(I) == I) NSGROUPS = NSGROUPS + 1
      end do
      allocate (SGMEMBERS(NMEM))
      NMEM = 0
      do I = 1, NP
         if (ISSGMEMBER(I)) then
            NMEM = NMEM + 1
            SGMEMBERS(NMEM) = I
         end if
      end do
      do I = 1, NP
         GPARENT(I) = I
      end do

      !...IBTYPE 4/24 pairs solved at the new time level (WEIR_SEMI_IMPLICIT) couple their two nodes
      !   in the GWCE matrix (IBTYPE 64 pairs are static group members and are linked already)
      NWP = 0
      if (WEIR_SEMI_IMPLICIT) then
         do I = 1, NVEL
            if (LBCODEI(I) /= 4 .and. LBCODEI(I) /= 24) cycle
            if (NBV(I) > 0 .and. IBCONN(I) > 0) NWP = NWP + 1
         end do
      end if

      !...GWCE neighbour table
      if (NSGROUPS == 0 .and. NWP == 0) then
         MNEIG = MNEI
         allocate (NeiTabG(MNP, MNEIG), NNeighG(MNP))
         NeiTabG = NeiTab
         NNeighG = NNeigh
         return
      end if

      allocate (CNT(NP), UNION(NP))
      allocate (TAB(NP, MNEI + 8))
      TAB = 0
      do I = 1, NP
         CNT(I) = NNeigh(I)
         TAB(I, 1:NNeigh(I)) = NeiTab(I, 1:NNeigh(I))
      end do
      do R = 1, NP
         if (FIRST(R) == 0) cycle
         !...UNION: the members and all their neighbours
         NU = 0
         I = FIRST(R)
         do while (I /= 0)
            call ADDUNIQUE(UNION, NU, I)
            do J = 1, NNeigh(I)
               call ADDUNIQUE(UNION, NU, NeiTab(I, J))
            end do
            I = NEXT(I)
         end do
         !...members of other groups next to this one (adjacent cross-sections of a condensed
         !   channel share elements): add their whole groups, since a row of this group can need
         !   the column of any of their members as representative
         NU0 = NU
         do K = 1, NU0
            N = UNION(K)
            if (ISSGMEMBER(N) .and. ROOT(N) /= R) then
               I = FIRST(ROOT(N))
               do while (I /= 0)
                  call ADDUNIQUE(UNION, NU, I)
                  I = NEXT(I)
               end do
            end if
         end do
         do K = 1, NU
            N = UNION(K)
            if (ISSGMEMBER(N) .and. ROOT(N) == R) then
               !...a member lists the whole union
               do L = 1, NU
                  call ADDTAB(N, UNION(L))
               end do
            else
               !...a neighbour lists all members
               I = FIRST(R)
               do while (I /= 0)
                  call ADDTAB(N, I)
                  I = NEXT(I)
               end do
            end if
         end do
      end do
      if (NWP > 0) then
         do I = 1, NVEL
            if (LBCODEI(I) /= 4 .and. LBCODEI(I) /= 24) cycle
            if (NBV(I) <= 0 .or. IBCONN(I) <= 0) cycle
            call ADDTAB(NBV(I), IBCONN(I))
            call ADDTAB(IBCONN(I), NBV(I))
         end do
      end if
      NCOL = maxval(CNT)
      MNEIG = max(MNEI, NCOL)
      allocate (NeiTabG(MNP, MNEIG), NNeighG(MNP))
      NeiTabG = 0
      NNeighG = 0
      do I = 1, NP
         NNeighG(I) = CNT(I)
         NeiTabG(I, 1:CNT(I)) = TAB(I, 1:CNT(I))
         !...unused slots point to the node itself (zero coefficient)
         NeiTabG(I, CNT(I) + 1:MNEIG) = I
      end do
      do IE = 1, NE
         SGELEM(IE) = ISSGMEMBER(NM(IE, 1)) .or. ISSGMEMBER(NM(IE, 2)) .or. ISSGMEMBER(NM(IE, 3))
      end do
      deallocate (TAB, CNT, UNION, ROOT, SIZEOF, FIRST, NEXT)

      write (MSG, '(a,i0,a,i0,a,i0,a,i0,a,i0,a)') "Consistent GWCE: ", NSGROUPS, &
         " static coupled groups (VEW pairs and condensed nodes) with ", NMEM, &
         " nodes, ", NWP, " IBTYPE 4/24 pair entries linked; GWCE neighbour table width ", MNEIG, &
         " (mesh table ", MNEI, ")."
      call allMessage(INFO, trim(MSG))

   contains

      !...append V to row N of TAB unless present; widen TAB when full
      subroutine ADDTAB(N, V)
         integer, intent(IN) :: N, V
         integer :: M
         integer, allocatable :: TMP(:, :)
         do M = 1, CNT(N)
            if (TAB(N, M) == V) return
         end do
         if (CNT(N) == size(TAB, 2)) then
            allocate (TMP(size(TAB, 1), 2*size(TAB, 2)))
            TMP = 0
            TMP(:, 1:size(TAB, 2)) = TAB
            call move_alloc(TMP, TAB)
         end if
         CNT(N) = CNT(N) + 1
         TAB(N, CNT(N)) = V
      end subroutine ADDTAB

      subroutine ADDUNIQUE(LIST, N, V)
         integer, intent(INOUT) :: LIST(:)
         integer, intent(INOUT) :: N
         integer, intent(IN) :: V
         integer :: M
         do M = 1, N
            if (LIST(M) == V) return
         end do
         N = N + 1
         LIST(N) = V
      end subroutine ADDUNIQUE

   end subroutine SETUP_COUPLED_GROUPS

   !-------------------------------------------------------------
   !  UPDATE_ACTIVE_GROUPS: active groups of this step and their representatives (GMAP). Two nodes of
   !  an IBTYPE=64 pair are joined if the pair is coupled (ISSUBMERGED64 /= 0) and both are wet, the
   !  wet members of a condensed group are joined; the representative is the member with the lowest
   !  node number. CHANGED is .true. if GMAP differs from the previous call's.
   !-------------------------------------------------------------
   subroutine UPDATE_ACTIVE_GROUPS(CHANGED)
      use GLOBAL, only: NODECODE
      use BOUNDARIES, only: NVEL, NBV, IBCONN, LBCODEI, ISSUBMERGED64
      use NodalAttributes, only: LoadCondensedNodes, NListCondensedNodes, &
                                 NNodesListCondensedNodes, ListCondensedNodes
      implicit none
      logical, intent(OUT) :: CHANGED
      integer :: I, K, L, M, N1, N2, R

      CHANGED = .false.
      if (NSGROUPS == 0) return
      do M = 1, size(SGMEMBERS)
         I = SGMEMBERS(M)
         GPARENT(I) = I
      end do
      do I = 1, NVEL
         if (LBCODEI(I) /= 64) cycle
         if (ISSUBMERGED64(I) == 0) cycle
         N1 = NBV(I)
         N2 = IBCONN(I)
         if (N1 <= 0 .or. N2 <= 0) cycle
         if (NODECODE(N1) == 0 .or. NODECODE(N2) == 0) cycle
         call GJOIN(N1, N2)
      end do
      if (LoadCondensedNodes) then
         do K = 1, NListCondensedNodes
            N1 = 0
            do L = 1, NNodesListCondensedNodes(K)
               N2 = ListCondensedNodes(K, L)
               if (N2 <= 0) cycle
               if (NODECODE(N2) == 0) cycle
               if (N1 == 0) then
                  N1 = N2
               else
                  call GJOIN(N1, N2)
               end if
            end do
         end do
      end if
      !...representative: the member with the lowest full-domain node number, so that every
      !   subdomain holding the group picks the same node (adcprep keeps groups whole)
      do M = 1, size(SGMEMBERS)
         I = SGMEMBERS(M)
         GCNT(GFIND(I)) = 0
      end do
      do M = 1, size(SGMEMBERS)
         I = SGMEMBERS(M)
         R = GFIND(I)
         if (GCNT(R) == 0) then
            GCNT(R) = I
         elseif (FULLNODE(I) < FULLNODE(GCNT(R))) then
            GCNT(R) = I
         end if
      end do
      do M = 1, size(SGMEMBERS)
         I = SGMEMBERS(M)
         GMAP(I) = GCNT(GFIND(I))
         if (GMAP(I) /= GMAPPREV(I)) CHANGED = .true.
         GMAPPREV(I) = GMAP(I)
      end do
   end subroutine UPDATE_ACTIVE_GROUPS

   !-------------------------------------------------------------
   !  GROUP_LOAD_SUM: sum the GWCE loads of each active group into its representative (zero at the
   !  other members) and replace ETA1 of the members by the group mean, as the lumped GWCE does.
   !  Sets the initial guess ETAS of the members other than the representative to zero.
   !-------------------------------------------------------------
   subroutine GROUP_LOAD_SUM(GWCE_LV, ETA1, ETAS)
      implicit none
      real(8), intent(INOUT) :: GWCE_LV(:), ETA1(:), ETAS(:)
      integer :: I, M, R

      if (NSGROUPS == 0) return
      do M = 1, size(SGMEMBERS)
         I = SGMEMBERS(M)
         GSUM(GMAP(I)) = 0.d0
         GCNT(GMAP(I)) = 0
      end do
      do M = 1, size(SGMEMBERS)
         I = SGMEMBERS(M)
         R = GMAP(I)
         GSUM(R) = GSUM(R) + ETA1(I)
         GCNT(R) = GCNT(R) + 1
         if (R /= I) then
            GWCE_LV(R) = GWCE_LV(R) + GWCE_LV(I)
            GWCE_LV(I) = 0.d0
            ETAS(I) = 0.d0
         end if
      end do
      do M = 1, size(SGMEMBERS)
         I = SGMEMBERS(M)
         R = GMAP(I)
         if (GCNT(R) > 1) ETA1(I) = GSUM(R)/dble(GCNT(R))
      end do
   end subroutine GROUP_LOAD_SUM

   !-------------------------------------------------------------
   !  GROUP_COPY_INCREMENT: after JCG, give the members of each active group the increment of its
   !  representative.
   !-------------------------------------------------------------
   subroutine GROUP_COPY_INCREMENT(ETAS)
      implicit none
      real(8), intent(INOUT) :: ETAS(:)
      integer :: I, M

      if (NSGROUPS == 0) return
      do M = 1, size(SGMEMBERS)
         I = SGMEMBERS(M)
         if (GMAP(I) /= I) ETAS(I) = ETAS(GMAP(I))
      end do
   end subroutine GROUP_COPY_INCREMENT

   !-------------------------------------------------------------
   !  GROUP_BLOCK_RESET: restore the GWCE matrix entries changed by GROUP_BLOCK_ADD. ASSEMBLED is
   !  .true. if the matrix was assembled again since (the saved values are then stale and dropped).
   !-------------------------------------------------------------
   subroutine GROUP_BLOCK_RESET(COEF, ASSEMBLED)
      implicit none
      real(8), intent(INOUT) :: COEF(:, :)
      logical, intent(IN) :: ASSEMBLED
      integer :: K
      if (.not. ASSEMBLED) then
         do K = NBLK, 1, -1
            COEF(BLKROW(K), BLKSLOT(K)) = BLKORIG(K)
         end do
      end if
      NBLK = 0
   end subroutine GROUP_BLOCK_RESET

   !-------------------------------------------------------------
   !  GROUP_BLOCK_ADD: add V to the GWCE matrix entry (row IROW, column ICOL), remembering the entry's
   !  value before its first change since the last GROUP_BLOCK_RESET
   !-------------------------------------------------------------
   subroutine GROUP_BLOCK_ADD(COEF, IROW, ICOL, V)
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      real(8), intent(INOUT) :: COEF(:, :)
      integer, intent(IN) :: IROW, ICOL
      real(8), intent(IN) :: V
      integer :: J, K
      integer, allocatable :: ITMP(:)
      real(8), allocatable :: RTMP(:)
      J = GROUP_SLOT(IROW, ICOL)
      if (J == 0) call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                                 message="GWCE: weir pair column missing from the GWCE neighbour table")
      do K = 1, NBLK
         if (BLKROW(K) == IROW .and. BLKSLOT(K) == J) exit
      end do
      if (K > NBLK) then
         if (.not. allocated(BLKROW)) then
            allocate (BLKROW(64), BLKSLOT(64), BLKORIG(64))
         elseif (NBLK == size(BLKROW)) then
            allocate (ITMP(2*NBLK)); ITMP(1:NBLK) = BLKROW(1:NBLK); call move_alloc(ITMP, BLKROW)
            allocate (ITMP(2*NBLK)); ITMP(1:NBLK) = BLKSLOT(1:NBLK); call move_alloc(ITMP, BLKSLOT)
            allocate (RTMP(2*NBLK)); RTMP(1:NBLK) = BLKORIG(1:NBLK); call move_alloc(RTMP, BLKORIG)
         end if
         NBLK = NBLK + 1
         BLKROW(NBLK) = IROW
         BLKSLOT(NBLK) = J
         BLKORIG(NBLK) = COEF(IROW, J)
      end if
      COEF(IROW, J) = COEF(IROW, J) + V
   end subroutine GROUP_BLOCK_ADD

   !-------------------------------------------------------------
   !  GROUP_SLOT: column of node ICOL in row IROW of the GWCE matrix (NeiTabG), 0 if absent
   !-------------------------------------------------------------
   integer function GROUP_SLOT(IROW, ICOL)
      implicit none
      integer, intent(IN) :: IROW, ICOL
      integer :: J
      GROUP_SLOT = 0
      if (IROW == ICOL) then
         GROUP_SLOT = 1
         return
      end if
      do J = 2, NNeighG(IROW)
         if (NeiTabG(IROW, J) == ICOL) then
            GROUP_SLOT = J
            return
         end if
      end do
   end function GROUP_SLOT

   !-------------------------------------------------------------
   !  FULLNODE: full-domain node number of subdomain node N (N in serial runs)
   !-------------------------------------------------------------
   integer function FULLNODE(N)
      use GLOBAL, only: NODES_LG
      implicit none
      integer, intent(IN) :: N
      FULLNODE = N
      if (allocated(NODES_LG)) FULLNODE = abs(NODES_LG(N))
   end function FULLNODE

   !-------------------------------------------------------------
   !  GROUP_IMAX: maximum of K over all ranks (K in serial runs)
   !-------------------------------------------------------------
   integer function GROUP_IMAX(K)
#ifdef CMPI
      use mpi_f08, only: MPI_Allreduce, MPI_INTEGER, MPI_MAX
      use GLOBAL, only: COMM
#endif
      implicit none
      integer, intent(IN) :: K
      GROUP_IMAX = K
#ifdef CMPI
      call MPI_Allreduce(K, GROUP_IMAX, 1, MPI_INTEGER, MPI_MAX, COMM)
#endif
   end function GROUP_IMAX

   !-------------------------------------------------------------
   !  union-find on GPARENT; the smaller node becomes the root
   !-------------------------------------------------------------
   integer function GFIND(N)
      implicit none
      integer, intent(IN) :: N
      integer :: R, K, NX
      R = N
      do while (GPARENT(R) /= R)
         R = GPARENT(R)
      end do
      K = N
      do while (GPARENT(K) /= R)
         NX = GPARENT(K)
         GPARENT(K) = R
         K = NX
      end do
      GFIND = R
   end function GFIND

   subroutine GJOIN(N1, N2)
      implicit none
      integer, intent(IN) :: N1, N2
      integer :: R1, R2
      R1 = GFIND(N1)
      R2 = GFIND(N2)
      if (R1 == R2) return
      if (R1 < R2) then
         GPARENT(R2) = R1
      else
         GPARENT(R1) = R2
      end if
   end subroutine GJOIN

end module vew1d
