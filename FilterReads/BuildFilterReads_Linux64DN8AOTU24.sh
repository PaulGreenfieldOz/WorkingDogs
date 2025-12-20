dotnet publish ./FilterReads.sln -c release /p:PublishProfile=Linux64DN8AOTFolderProfile.pubxml
mkdir -p Linux64DN8AOTU24
cp FilterReads/bin/Release/net8.0/publish/linux-x64/FilterReads ./Linux64DN8AOTU24/
