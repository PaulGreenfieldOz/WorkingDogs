dotnet publish ./CondenseProkkaTbl.sln -c release /p:PublishProfile=Linux64DN8AOTFolderProfile.pubxml
mkdir -p Linux64DN8AOTU24
cp CondenseProkkaTbl/bin/Release/net8.0/publish/linux-x64/CondenseProkkaTbl ./Linux64DN8AOTU24/
