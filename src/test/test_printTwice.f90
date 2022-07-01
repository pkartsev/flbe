subroutine printTwice_openFile( filename, message )
character(len=*) :: filename
character(len=*), optional :: message

	open( debugFileNumber, file=filename )
	if (present(message)) then
		call printTwice( message, .TRUE. )
	end if

end subroutine printTwice_openFile


subroutine printTwice_closeFile( message )
character(len=*), optional :: message

	if (present(message)) then
!		call printTwice( message//CHAR(10), .TRUE. )
		call printTwice( message, .TRUE. )
	end if
	if (verbose) then
		close( debugFileNumber )
	else
		close( debugFileNumber, status='delete' )
	end if
	
end subroutine printTwice_closeFile


subroutine printTwice( s, forcePrint )
character(len=*) :: s
logical, optional :: forcePrint
logical doPrint
	
	doPrint = verbose
	if (present(forcePrint)) doPrint = forcePrint
	
	write( debugFileNumber, '(A)') s
!	if (doPrint) write(*, '(A)') s
	if (doPrint) then
		write(*, '(A)') s
	end if
	
end subroutine printTwice


subroutine printTwiceSL( s, L, forcePrint )
character(len=*) :: s
logical L
logical, optional :: forcePrint

logical doPrint

	doPrint = verbose
	if (present(forcePrint)) doPrint = forcePrint
	
	write( debugFileNumber, '(A,L)') S, L
	if (doPrint) write(*, '(A,L)') s, L
	
end subroutine printTwiceSL


subroutine printTwiceSIS( s1, i, s2, forcePrint )
integer i
character(len=*) :: s1, s2
logical, optional :: forcePrint

logical doPrint
	
	doPrint = verbose
	if (present(forcePrint)) doPrint = forcePrint
	
	write( debugFileNumber, '(A,I0,A)') s1, i, s2
	if (doPrint) write(*, '(A,I0,A)') s1, i, s2
	
end subroutine printTwiceSIS


subroutine printTwiceSFS( s1, f, s2 )
integer i
double precision f
character(len=*) :: s1, s2

	write( debugFileNumber, '(A,ES18.8,A)') s1,f,s2
	if (verbose) write(*, '(A,ES18.8,A)') s1,f,s2

end subroutine printTwiceSFS


subroutine printTwiceSFSFS( s1, f1, s2, f2, s3 )
integer i
double precision f1, f2
character(len=*) :: s1, s2, s3

	write( debugFileNumber, '(A,ES18.8,A,ES18.8,A)') s1,f1,s2,f2,s3
	if (verbose) write(*, '(A,ES18.8,A,ES18.8,A)') s1,f1,s2,f2,s3

end subroutine printTwiceSFSFS

