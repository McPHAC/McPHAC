


      subroutine opalrho(p,t6,r)
      
c     This is the routine called by the C program, to interpolate
c     values from the OPAL EOS file;  Temperature (in units of 1e6 K)
c     and pressure P (units???) are input, and R=rho/T6**3 is returned
      
      parameter (mx=6,mv=10,nr=169,nt=191)
      character*3 fixedTP
      data fixedTP/'yes'/       ! 'yes' for fixed T6,P; 'no' for fixed T6,rho
      common/eeos/esact,eos(mv)
      iorder=1                  ! only interested in density; see esac
      irad=0                    ! does not add radiation corrections
      
      x=1.0                     ! pure Hydrogen

c      write (*,'("T6, p, x=",3e14.4)') t6,p,x 
      if (fixedTP .eq. 'yes') then
         r=rhoofp (x,t6,p,irad) ! calculate density (r) for fixed P. 
      endif 
c     call esac(x,t6,r,iorder,irad) ! calc EOS; use r from rhoop
c     don't need the above statment if all I need is density! 
      return
      end
      
