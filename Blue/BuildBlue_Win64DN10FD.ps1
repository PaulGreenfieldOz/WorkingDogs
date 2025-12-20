dotnet publish .\Blue.sln -c release /p:PublishProfile=Win64DN10FDFolderProfile.pubxml
mkdir -force Win64DN10FD
copy Blue\bin\Release\net10.0\publish\win-x64\Blue.exe .\Win64DN10FD\
