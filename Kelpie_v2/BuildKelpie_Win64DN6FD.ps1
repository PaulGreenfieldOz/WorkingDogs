dotnet publish .\Kelpie_v2.csproj -c release /p:PublishProfile=Win64DN6FDFolderProfile.pubxml
mkdir Win64DN6FD
copy bin\Release\net6.0\publish\win-x64\Kelpie_v2.exe .\Win64DN6FD\
