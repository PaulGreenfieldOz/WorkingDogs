dotnet publish ./Tessel.csproj -c release /p:PublishProfile=Linux64DN6FDFolderProfile.pubxml
mkdir -p Linux64DN6FD
cp bin/Release/net6.0/publish/linux-x64/Tessel ./Linux64DN6FD/
