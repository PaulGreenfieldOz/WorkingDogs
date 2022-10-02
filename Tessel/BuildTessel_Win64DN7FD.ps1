dotnet publish .\Tessel.csproj -c release /p:PublishProfile=Win64DN7FDFolderProfile.pubxml
mkdir Win64DN7FD
copy bin\Release\net7.0\publish\win-x64\Tessel.exe .\Win64DN7FD\
