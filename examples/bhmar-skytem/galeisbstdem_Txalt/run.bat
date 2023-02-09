@echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%
set path=..\..\..\third_party\gdal;%path%
set path=..\..\..\third_party\vs2015;%path%

REM Standalone
REM galeisbstdem.exe galeisbstdem.con

REM If you don't have MPI installed use the non-MPI version
REM Use n OpenMP threads 
galeisbstdem-nompi.exe galeisbstdem.con 3

pause


