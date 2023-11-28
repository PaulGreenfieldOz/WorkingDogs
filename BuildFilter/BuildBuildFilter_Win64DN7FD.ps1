dotnet publish .\BuildFilter.csproj -c release /p:PublishProfile=Win64DN7FDFolderProfile.pubxml
mkdir -force Win64DN7FD
copy bin\Release\net7.0\publish\win-x64\BuildFilter.exe .\Win64DN7FD\
