dotnet publish ./Clancy.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml
mkdir -p Linux64DN10AOTU24
cp Clancy/bin/Release/net10.0/publish/linux-x64/Clancy ./Linux64DN10AOTU24/
