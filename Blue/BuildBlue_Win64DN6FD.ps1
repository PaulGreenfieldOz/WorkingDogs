dotnet publish .\Blue.csproj -c release /p:PublishProfile=Win64DN6FDFolderProfile.pubxml
mkdir -force Win64DN6FD
copy bin\Release\net6.0\publish\win-x64\Blue.exe .\Win64DN6FD\
