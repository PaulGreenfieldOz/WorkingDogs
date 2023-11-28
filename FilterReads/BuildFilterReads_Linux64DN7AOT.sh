dotnet publish ./FilterReads.csproj -c release /p:PublishProfile=Linux64DN7AOTFolderProfile.pubxml
mkdir -p Linux64DN7AOT
cp bin/Release/net7.0/publish/linux-x64/FilterReads ./Linux64DN7AOT/
