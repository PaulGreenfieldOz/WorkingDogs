dotnet publish ./BuildFilter.csproj -c release /p:PublishProfile=Linux64DN6FDFolderProfile.pubxml
mkdir Linux64DN6FD
cp bin/Release/net6.0/publish/linux-x64/BuildFilter ./Linux64DN6FD/
