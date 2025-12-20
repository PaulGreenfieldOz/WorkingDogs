dotnet publish ./Blue.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml
mkdir -p Linux64DN10AOTU24
cp Blue/bin/Release/net10.0/publish/linux-x64/Blue ./Linux64DN10AOTU24/
