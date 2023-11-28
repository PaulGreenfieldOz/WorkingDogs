dotnet publish .\BuildFilter.csproj -c release /p:PublishProfile=Win64DN8FDFolderProfile.pubxml
mkdir -force Win64DN8FD
copy bin\Release\net8.0\publish\win-x64\BuildFilter.exe .\Win64DN8FD\
