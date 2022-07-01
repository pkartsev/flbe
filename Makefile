
DEVELOP=no
release_name := 101-3-rc

#------------------------------------------------------------------------------------
ifeq ($(DEVELOP),yes)
	COMPILATION := -fbacktrace -fcheck=all -g3 -fopenmp 
else
	COMPILATION := -O3 -march=native -mavx -funsafe-math-optimizations -fopenmp
endif

# https://gcc.gnu.org/onlinedocs/gfortran/Error-and-Warning-Options.html
#
# -Wall 
# This currently includes -Waliasing, -Wampersand, -Wconversion, -Wsurprising, 
# -Wc-binding-type, -Wintrinsics-std, -Wtabs, -Wintrinsic-shadow, -Wline-truncation,
# -Wtarget-lifetime, -Winteger-division, -Wreal-q-constant, -Wunused and -Wundefined-do-loop
#

ifeq ($(OS),Windows_NT)
    MAIN_EXE := main.exe
    DEL := del
	FFT_PATH := -I D:/development/fftw/ -L D:/development/fftw/ -lfftw3-3
	MOD_PATH := -J $(TEMP)
else
    MAIN_EXE := ./main.exe
    DEL := rm -f
	FFT_PATH := -I /usr/include -lfftw3 -lfftw3_omp
	MOD_PATH := -J /tmp/
endif

FLBE_PATH := -Isrc/flbe 
WORK_FILES := *.exe *.tex *.out *.log *.aux *.txt
#------------------------------------------------------------------------------------

OPTS:= $(FFT_PATH) $(COMPILATION) $(MOD_PATH) $(FLBE_PATH)

#------------------------------------------------------------------------------------

all: minimal.exe
#all: tests.exe
#all: benchmark_v0.exe benchmark.exe
#all: main.exe
#all: tests.exe benchmark.exe benchmark_v0.exe main.exe

1.pdf: 1.tex
		pdflatex 1.tex

1.tex: main.exe
		$(MAIN_EXE)

main.exe: src/demo/main.f90 src/flbe/*.f90 Makefile
		gfortran src/demo/main.f90 $(OPTS) -o main.exe

tests.exe: src/test/*.f90 src/flbe/*.f90 Makefile
		gfortran src/test/tests.f90 $(OPTS) -Isrc/test -o tests.exe

benchmark.exe: src/benchmark/benchmark.f90 src/flbe/*.f90 Makefile
		gfortran src/benchmark/benchmark.f90 $(OPTS) -o benchmark.exe

benchmark_v0.exe: src/test/benchmark_v0.f90 src/flbe/*.f90 src/test/v0*.f90 Makefile
		gfortran src/test/benchmark_v0.f90 $(OPTS) -o benchmark_v0.exe

minimal.exe: src/demo/minimal.f90 src/flbe/*.f90 Makefile
		gfortran src/demo/minimal.f90 $(OPTS) -o minimal.exe

clean:
		$(DEL) $(WORK_FILES) *.pdf *.html *.dat $(release_name).zip

zip:
		zip $(release_name).zip Makefile src/*/*.f90 *.txt *.png *.pdf
