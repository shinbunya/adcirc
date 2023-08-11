      MODULE MP_INTERFACE
      use sizes, only: SZ
      use global, only: DEBUG, INFO, ECHO, WARNING, ERROR,
     &   scratchMessage, logMessage, screenMessage, allMessage,
     &   setMessageSource, openFileForRead, unsetMessageSource,
     &   scratchFormat
      
      IMPLICIT NONE
      
      SAVE

      CONTAINS

C======================================================================
      SUBROUTINE MPINTF_Pack_EdgeValue()
C----------------------------------------------------------------------
      USE GLOBAL, ONLY: ETA1, QX1, QY1
      USE MESH, ONLY: DP
      USE BOUNDARIES, ONLY: NVEL, NBV, LBCODEI
      USE MP_MESSENGER, ONLY: MPMSG_BUFF_VAL, nEdge
      IMPLICIT NONE

      INTEGER :: IEDG, J1, J2, NBDJ1, NBDJ2
      REAL(SZ) :: ZE, HB, QX, QY
      CHARACTER(2024) :: message

      call setMessageSource("MP_Pack_EdgeValue")

      IEDG = 0
      MPMSG_BUFF_VAL(:) = 0.D0
      
      DO J1=1,NVEL-1
         J2 = J1 + 1
         NBDJ1 = NBV(J1)
         NBDJ2 = NBV(J2)
         IF(LBCODEI(J1).EQ.22.AND.LBCODEI(J2).EQ.22) THEN
            ZE = 0.5d0*(ETA1(NBDJ1) + ETA1(NBDJ2))
            HB = 0.5d0*(DP(NBDJ1) + DP(NBDJ2))
            QX = 0.5d0*(QX1(NBDJ1) + QX1(NBDJ2))
            QY = 0.5d0*(QY1(NBDJ1) + QY1(NBDJ2))
            MPMSG_BUFF_VAL(IEDG*4 + 1) = ZE
            MPMSG_BUFF_VAL(IEDG*4 + 2) = HB
            MPMSG_BUFF_VAL(IEDG*4 + 3) = QX
            MPMSG_BUFF_VAL(IEDG*4 + 4) = QY
            IEDG = IEDG + 1
         ENDIF
      END DO

      IF(IEDG.NE.nEdge) THEN
         WRITE(message,'(A,I8,A,I8)') ' - We got ', IEDG, ' while expecting ', nEdge
         call allMessage(ERROR, 'Inconsistent number of edges')
         call allMessage(ERROR, trim(message))
         call terminate()
      ENDIF

      call unsetMessageSource()

      RETURN
      END SUBROUTINE MPINTF_Pack_EdgeValue
      

C======================================================================
      SUBROUTINE MPINTF_Unpack_FluxValue()
C----------------------------------------------------------------------
      IMPLICIT NONE
      ! Receive flux values from Coupler and unpack them to QN1
      RETURN
      END SUBROUTINE MPINTF_Unpack_FluxValue


C======================================================================
      SUBROUTINE terminate(NO_MPI_FINALIZE)
C----------------------------------------------------------------------
#if defined(CMPI) || defined(CMPMD)
      USE MESSENGER
#endif
      IMPLICIT NONE
C
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
C
      call setMessageSource("terminate")
#if defined(MP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#if defined(CMPI) || defined(CMPMD)
      subdomainFatalError = .true.
      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
      ELSE
        CALL MSG_FINI()
      ENDIF
#endif
      CALL EXIT(1)
C
#if defined(MP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()

      END SUBROUTINE terminate      
      END MODULE MP_INTERFACE