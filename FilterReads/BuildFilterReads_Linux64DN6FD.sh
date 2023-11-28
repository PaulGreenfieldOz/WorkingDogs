dotnet publish ./FilterReads.csproj -c release /p:PublishProfile=Linux64DN6FDFolderProfile.pubxml
mkdir -p Linux64DN6FD
cp bin/Release/net6.0/publish/linux-x64/FilterReads ./Linux64DN6FD/
