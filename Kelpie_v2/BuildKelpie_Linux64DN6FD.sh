dotnet publish ./Kelpie_v2.csproj -c release /p:PublishProfile=Linux64DN6FDFolderProfile.pubxml
mkdir Linux64DN6FD
cp bin/Release/net6.0/publish/linux-x64/Kelpie_v2 ./Linux64DN6FD/
