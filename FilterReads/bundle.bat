REM https://download.mono-project.com/runtimes/raw/

mkdir OSX-10.9
mkdir ubuntu-16.04
mkdir ubuntu-18.04
mkdir debian-8
mkdir debian-9
mkdir Windows

call mkbundle -o OSX-10.9\FilterReads FilterReads.exe --simple --cross mono-5.20.1-osx-10.9-x64
call mkbundle -o ubuntu-16.04\FilterReads FilterReads.exe --simple --cross mono-5.20.1-ubuntu-16.04-x64
call mkbundle -o ubuntu-18.04\FilterReads FilterReads.exe --simple --cross mono-5.20.1-ubuntu-18.04-x64
call mkbundle -o debian-8\FilterReads FilterReads.exe --simple --cross mono-5.20.1-debian-8-x64
call mkbundle -o debian-9\FilterReads FilterReads.exe --simple --cross mono-5.20.1-debian-9-x64
copy FilterReads.exe Windows
