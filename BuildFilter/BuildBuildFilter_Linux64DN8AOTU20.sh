dotnet publish ./BuildFilter.sln -c release /p:PublishProfile=Linux64DN8AOTFolderProfile.pubxml
mkdir -p Linux64DN8AOTU20
cp BuildFilter/bin/Release/net8.0/publish/linux-x64/BuildFilter ./Linux64DN8AOTU20/
