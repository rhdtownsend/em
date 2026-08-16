      subroutine mpi_fitness (num_jobs, npar, oldph, fitness)
c ---------------------------------------------
c     parallel fitness evaluation using MPI
c ---------------------------------------------
      implicit none

      include 'mpif.h'
      
      integer ndone, nproc, ierr, nspawn, valuest(8)
      integer msgtype, job, num_jobs, trial, slave
      integer par, npar, status(MPI_STATUS_SIZE)
      double precision data(32),result,age
      double precision R_Ro,Rcz,Teff,M_H,Dcz
      real oldph(32,1024), fitness(1024), youth(1024)
      character*10 datest, timest, zonest
      logical receiving

c ---------------------------------------------
c     initialize counter
c ---------------------------------------------
      ndone = 0

c ---------------------------------------------
c     determine number of processors
c ---------------------------------------------
      call mpi_comm_size( MPI_COMM_WORLD, nproc, ierr )

c ---------------------------------------------
c     run jobs on slave nodes only
c ---------------------------------------------
      nspawn = nproc-1

c ---------------------------------------------
c     send an initial job to each node
c ---------------------------------------------
      do job=1,nspawn
         trial = job
         slave = job
         call sendjob(trial,slave,npar,oldph)
    8    format("job",i3.3," ->",4(1x,f4.2))
	 write(*,8) trial,oldph(1,trial),oldph(2,trial),
     +              oldph(3,trial),oldph(4,trial)
      enddo
      write(*,*)
      call date_and_time(datest, timest, zonest, valuest)
      write(*,'(A10,": started job000")') timest
      timest = 'PARAMETERS'

      do job=1,num_jobs
c ---------------------------------------------
c     listen for responses
c ---------------------------------------------
 25      msgtype  = 2 
         call mpi_iprobe( MPI_ANY_SOURCE, msgtype, MPI_COMM_WORLD,
     +                    receiving, status, ierr )

         if (receiving) then
c ---------------------------------------------
c     get data from responding node
c ---------------------------------------------
            call mpi_recv( slave, 1, MPI_INTEGER, MPI_ANY_SOURCE,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( trial, 1, MPI_INTEGER, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( result, 1, MPI_DOUBLE_PRECISION, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( age, 1, MPI_DOUBLE_PRECISION, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( R_Ro, 1, MPI_DOUBLE_PRECISION, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( Rcz, 1, MPI_DOUBLE_PRECISION, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( Teff, 1, MPI_DOUBLE_PRECISION, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( M_H, 1, MPI_DOUBLE_PRECISION, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
c            call mpi_recv( fn0, 1, MPI_DOUBLE_PRECISION, slave,
c     +                     msgtype, MPI_COMM_WORLD, status, ierr )
c            call mpi_recv( Dnu_calc, 1, MPI_DOUBLE_PRECISION, slave,
c     +                     msgtype, MPI_COMM_WORLD, status, ierr )
c            call mpi_recv( chisq, 4, MPI_DOUBLE_PRECISION, slave,
c     +                     msgtype, MPI_COMM_WORLD, status, ierr )

            Dcz = Rcz/R_Ro
            youth(trial) = age
            fitness(trial) = result
            ndone = ndone + 1

c remove system call to improve performance
            call date_and_time(datest, timest, zonest, valuest)
c    9       format(A10,": job",i3.3,4(1x,f4.2),2(1x,e12.6),1x,
c     +        F5.3,1x,F7.5,1x,F7.1,1x,F6.3,1x,F7.2,1x,F7.2,4(1x,F7.3))
c	    write(*,9) timest,trial,oldph(1,trial),oldph(2,trial),
c     +        oldph(3,trial),oldph(4,trial),fitness(trial),age,
c     +        R_Ro,Dcz,Teff,M_H,fn0,Dnu_calc,chisq(1),chisq(2),
c     +        chisq(3),chisq(4)
    9       format(A10,": job",i3.3,4(1x,f4.2),2(1x,e12.6),1x,
     +         F5.3,1x,F7.5,1x,F7.1,1x,F6.3)
	    write(*,9) timest,trial,oldph(1,trial),oldph(2,trial),
     +        oldph(3,trial),oldph(4,trial),fitness(trial),age,
     +        R_Ro,Dcz,Teff,M_H

c ---------------------------------------------
c     send new job to responding node
c ---------------------------------------------
 140        if (ndone .LE. (num_jobs-nspawn)) then
               trial = job + nspawn
               call sendjob(trial,slave,npar,oldph)
	       write(*,8) trial,oldph(1,trial),oldph(2,trial),
     +                    oldph(3,trial),oldph(4,trial)
            endif
            goto 100
         endif

c ---------------------------------------------
c     return to listen again or move on
c ---------------------------------------------
         if (.NOT. receiving) goto 25

         goto 199
 100     continue
      enddo

c ---------------------------------------------
c     ready for next generation of jobs
c ---------------------------------------------

c      open(55,file='param.log',status='unknown')
c      do job=1,num_jobs
c 198     format(i3.3,4(1x,f4.2),2(1x,e12.6))
c         write(55,198) job,oldph(1,job),oldph(2,job),oldph(3,job),
c     +        oldph(4,job),youth(job),fitness(job)
c      enddo
c      close(55)

 199  return
      end

c**********************************************************************
      subroutine sendjob(trial,slave,npar,oldph)

      implicit none

      include 'mpif.h'

      integer trial, slave, par, npar, ierr, msgtype
      double precision data(32)
      real oldph(32,1024)

      do par=1,npar
         data(par) = oldph(par,trial)
      enddo

      msgtype = 1
      call mpi_send( trial, 1, MPI_INTEGER, slave,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( npar, 1, MPI_INTEGER, slave,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( data, npar, MPI_DOUBLE_PRECISION, slave,
     +               msgtype, MPI_COMM_WORLD, ierr )

      return
      end

c**********************************************************************
