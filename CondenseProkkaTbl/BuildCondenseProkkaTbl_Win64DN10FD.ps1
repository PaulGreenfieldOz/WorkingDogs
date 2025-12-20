dotnet publish .\CondenseProkkaTbl.sln -c release /p:PublishProfile=Win64DN10FDFolderProfile.pubxml
mkdir -force Win64DN10FD
copy CondenseProkkaTbl\bin\Release\net10.0\publish\win-x64\CondenseProkkaTbl.exe .\Win64DN10FD\
