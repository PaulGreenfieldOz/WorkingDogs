dotnet publish .\FilterReads.csproj -c release /p:PublishProfile=Win64DN8FDFolderProfile.pubxml
mkdir Win64DN8FD
copy bin\Release\net8.0\publish\win-x64\FilterReads.exe .\Win64DN8FD\
