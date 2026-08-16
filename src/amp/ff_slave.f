      subroutine ff_slave 
c ---------------------------------------------
c fitness function slave program
c ---------------------------------------------
      use em_lib, only: age_m, R_m, Teff_m, FeH_m, Rcz_m
      implicit none

      include 'mpif.h'

      integer myid, ierr, status(MPI_STATUS_SIZE)
      integer master, msgtype, trial, n, i

      double precision data(32), result
      real userff
      external userff

c ---------------------------------------------
c identify this slave task
c ---------------------------------------------
      call mpi_comm_rank( MPI_COMM_WORLD, myid, ierr )

c ---------------------------------------------
c listen for a new job
c ---------------------------------------------
      master = 0
 25   msgtype = 1

c ---------------------------------------------
c receive data from master host
c ---------------------------------------------
      call mpi_recv( trial, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )
      if (trial .EQ. -1) goto 99
      call mpi_recv( n, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( data, n, MPI_DOUBLE_PRECISION, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )

      do i=1,n
        data(i) = INT((100.*data(i))+0.5)/100.d0
      enddo

c ---------------------------------------------
c perform calculations with data
c ---------------------------------------------
      result = userff( n, data, myid )

c ---------------------------------------------
c send result to master host
c ---------------------------------------------      
      msgtype = 2 
      call mpi_send( myid, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( trial, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( result, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( age_m, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( R_m, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( Rcz_m, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( Teff_m, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( FeH_m, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )

c ---------------------------------------------
c go back for more work
c ---------------------------------------------
      goto 25

 99   return
      end

c*********************************************************************
