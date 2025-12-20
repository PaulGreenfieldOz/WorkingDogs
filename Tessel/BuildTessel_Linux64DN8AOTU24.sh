dotnet publish ./Tessel.sln -c release /p:PublishProfile=Linux64DN8AOTFolderProfile.pubxml
mkdir -p Linux64DN8AOTU24
cp Tessel/bin/Release/net8.0/publish/linux-x64/Tessel ./Linux64DN8AOTU24/
