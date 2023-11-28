dotnet publish ./BuildFilter.csproj -c release /p:PublishProfile=Linux64DN8AOTFolderProfile.pubxml
mkdir -p Linux64DN8AOT
cp bin/Release/net8.0/publish/linux-x64/BuildFilter ./Linux64DN8AOT/
