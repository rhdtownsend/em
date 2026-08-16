      program mpi_pikaia
c ----------------------------------------------
c Front end for parallel PIKAIA to work with MPI
c ----------------------------------------------
      implicit none

      include 'mpif.h'

      integer ierr,myid,nproc,nslaves,rc
      integer trial,slave,msgtype

c ----------------------------------------------
c     Initialize MPI
c ----------------------------------------------
      call mpi_init( ierr )
      call mpi_comm_rank( MPI_COMM_WORLD, myid, ierr )
      call mpi_comm_size( MPI_COMM_WORLD, nproc, ierr )
      nslaves=nproc-1

c ----------------------------------------------
c     Master program (parallel PIKAIA)
c ----------------------------------------------
      if (myid .EQ. 0) then
         call pikaia_master
c ----------------------------------------------
c     Finish with shutdown signal to slaves
c ----------------------------------------------
         trial = -1
         msgtype = 1
         do slave=1,nslaves
            call mpi_send( trial, 1, MPI_INTEGER, slave,
     +                     msgtype, MPI_COMM_WORLD, ierr )
         enddo
c ----------------------------------------------
c     Slave tasks (fitness evaluation)
c ----------------------------------------------
      elseif (myid .GT. 0) then
         call ff_slave
      endif

      call mpi_finalize(rc)
      stop
      end
