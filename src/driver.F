C******************************************************************************
C******************************************************************************
C*    Non-ESMF ADCIRC DRIVER
C******************************************************************************
C******************************************************************************
      PROGRAM ADCIRC
C
      USE ADCIRC_Mod, ONLY : ADCIRC_Init, ADCIRC_Run, ADCIRC_Final
#ifdef CMPMD
      USE MPI
      USE MP_MESSENGER, ONLY : MPMSG_Init, MPMSG_Fini,
     &    COLOR_PG_CG, MPMSG_Split_Comm, COLOR_CPLD_IN, COLOR_CPLD_OUT,
     &    MPMSG_Set_EdgeBuffer
      USE BOUNDARIES, ONLY : nMPEdge
#endif
C
#ifdef CSWAN
Casey 090302: Include the following routines for coupling to unstructured SWAN.
      USE Couple2swan, ONLY: PADCSWAN_INIT, PADCSWAN_FINAL
      USE SIZES, ONLY: MNPROC,MYPROC
#elif defined(CMPMD)
      USE SIZES, ONLY: MNPROC,MYPROC
#endif
C
      IMPLICIT NONE
C
#ifdef CMPMD
      INTEGER :: WRANK, ADC_COMM, IERR
#endif
C
#ifdef CMPMD
C.......Initialize MPI and split COMM_WORLD for multi-program runs
      CALL MPMSG_Init(COLOR_PG_CG, ADC_COMM)

C.......Initialize ADCIRC nodes with ADCIRC communicator
      CALL ADCIRC_Init (ADC_COMM)

C.......Initialize MPI communicator among coupled nodes (= compute nodes + coupler nodes)
      IF(MYPROC.LT.MNPROC) THEN
         CALL MPMSG_Split_Comm(COLOR_CPLD_IN)
      ELSE
         CALL MPMSG_Split_Comm(COLOR_CPLD_OUT)
      ENDIF

C.......Set send/recv buffers for edge values
      CALL MPMSG_Set_EdgeBuffer(nMPEdge)

#else
C.......Initialize ADCIRC for single program with single or multiple nodes
      CALL ADCIRC_Init
#endif
C
#ifdef CSWAN
Casey 090302: Allow SWAN to initialize stuff before the start
C             of the time step loop.  This subroutine is inside
C             the 'couple2swan.F' src file.
Casey 110518: Added this IF statement.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_INIT
      ENDIF
#endif
C
      CALL ADCIRC_Run
C
#ifdef CSWAN
Casey 090302: Let SWAN clean up stuff.
Casey 110518: Added this IF statement.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_FINAL
      ENDIF
#endif
C
      CALL ADCIRC_Final
C
#ifdef CMPMD
C.......Finalize MPI
      CALL MPI_FINALIZE (IERR)
#endif
C
      CALL EXIT(0)
      END PROGRAM ADCIRC
