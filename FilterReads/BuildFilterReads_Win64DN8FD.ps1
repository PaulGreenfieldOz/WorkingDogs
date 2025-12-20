dotnet publish .\FilterReads.sln -c release /p:PublishProfile=Win64DN8FDFolderProfile.pubxml
mkdir -force Win64DN8FD
copy FilterReads\bin\Release\net8.0\publish\win-x64\FilterReads.exe .\Win64DN8FD\
