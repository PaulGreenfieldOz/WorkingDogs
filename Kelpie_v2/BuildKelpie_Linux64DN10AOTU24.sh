dotnet publish ./Kelpie_v2.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml
mkdir -p Linux64DN10AOTU24
cp Kelpie_v2/bin/Release/net10.0/publish/linux-x64/Kelpie_v2 ./Linux64DN10AOTU24/
cp Kelpie_v2/bin/Release/net10.0/publish/linux-x64/Kelpie_v2 ./Linux64DN10AOTU24/Kelpie
