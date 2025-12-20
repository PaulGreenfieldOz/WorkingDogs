dotnet publish .\CondenseProkkaTbl.sln -c release /p:PublishProfile=Win64DN8FDFolderProfile.pubxml
mkdir -force Win64DN8FD
copy CondenseProkkaTbl\bin\Release\net8.0\publish\win-x64\CondenseProkkaTbl.exe .\Win64DN8FD\
