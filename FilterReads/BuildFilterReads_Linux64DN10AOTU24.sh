dotnet publish ./FilterReads.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml
mkdir -p Linux64DN10AOTU24
cp FilterReads/bin/Release/net10.0/publish/linux-x64/FilterReads ./Linux64DN10AOTU24/
