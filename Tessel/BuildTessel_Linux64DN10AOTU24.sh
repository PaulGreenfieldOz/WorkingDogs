dotnet publish ./Tessel.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml
mkdir -p Linux64DN10AOTU24
cp Tessel/bin/Release/net10.0/publish/linux-x64/Tessel ./Linux64DN10AOTU24/
