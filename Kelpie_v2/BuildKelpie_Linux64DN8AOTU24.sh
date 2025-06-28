dotnet publish ./Kelpie_v2.csproj -c release /p:PublishProfile=Linux64DN8AOTFolderProfile.pubxml
mkdir -p Linux64DN8AOTU24
cp bin/Release/net8.0/publish/linux-x64/Kelpie_v2 ./Linux64DN8AOTU24/
