REM https://download.mono-project.com/runtimes/raw/

mkdir OSX-10.7
mkdir ubuntu-16.04
mkdir ubuntu-18.04
mkdir debian-8
mkdir debian-9

call mkbundle -o OSX-10.7\BuildFilter BuildFilter.exe --simple --cross mono-5.12.0-osx-10.7-x64
call mkbundle -o ubuntu-16.04\BuildFilter BuildFilter.exe --simple --cross mono-5.12.0-ubuntu-16.04-x64
call mkbundle -o ubuntu-18.04\BuildFilter BuildFilter.exe --simple --cross mono-5.12.0-ubuntu-18.04-x64
call mkbundle -o debian-8\BuildFilter BuildFilter.exe --simple --cross mono-5.12.0-debian-8-x64
call mkbundle -o debian-9\BuildFilter BuildFilter.exe --simple --cross mono-5.12.0-debian-9-x64

