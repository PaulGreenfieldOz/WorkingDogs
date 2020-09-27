REM https://download.mono-project.com/runtimes/raw/

mkdir OSX-10.7
mkdir ubuntu-16.04
mkdir ubuntu-18.04
mkdir debian-8
mkdir debian-9
mkdir Windows

call mkbundle -o OSX-10.7\Blue Blue.exe --simple --cross mono-5.14.0-osx-10.7-x64
call mkbundle -o ubuntu-16.04\Blue Blue.exe --simple --cross mono-5.14.0-ubuntu-16.04-x64
call mkbundle -o ubuntu-18.04\Blue Blue.exe --simple --cross mono-5.14.0-ubuntu-18.04-x64
call mkbundle -o debian-8\Blue Blue.exe --simple --cross mono-5.14.0-debian-8-x64
call mkbundle -o debian-9\Blue Blue.exe --simple --cross mono-5.14.0-debian-9-x64

copy Blue.exe Windows
