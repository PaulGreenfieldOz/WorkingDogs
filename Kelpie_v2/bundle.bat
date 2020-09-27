REM https://download.mono-project.com/runtimes/raw/

rem start mono command prompt and cd to this directory - and then run this batch file
rem new runtimes can be downloaded from above link - but need renaming to .zip files before running

mkdir OSX-10.9
mkdir ubuntu-16.04
mkdir ubuntu-18.04
mkdir debian-9
REM mkdir debian-10
mkdir Windows

call mkbundle -o OSX-10.9\Kelpie_v2 Kelpie_v2.exe --simple --cross mono-5.20.1-osx-10.9-x64
call mkbundle -o ubuntu-16.04\Kelpie_v2 Kelpie_v2.exe --simple --cross mono-5.20.1-ubuntu-16.04-x64
call mkbundle -o ubuntu-18.04\Kelpie_v2 Kelpie_v2.exe --simple --cross mono-5.20.1-ubuntu-18.04-x64
call mkbundle -o debian-9\Kelpie_v2 Kelpie_v2.exe --simple --cross mono-5.20.1-debian-9-x64
rem call mkbundle -o debian-10\Kelpie_v2 Kelpie_v2.exe --simple --cross mono-6.8.0-debian-10-x64
copy Kelpie_v2.exe Windows

