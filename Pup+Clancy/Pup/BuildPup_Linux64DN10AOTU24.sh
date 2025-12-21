dotnet publish ./Pup.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml
mkdir -p Linux64DN10AOTU24
cp Pup/bin/Release/net10.0/publish/linux-x64/Pup ./Linux64DN10AOTU24/
