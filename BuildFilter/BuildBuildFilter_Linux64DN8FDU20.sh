dotnet publish ./BuildFilter.sln -c release /p:PublishProfile=Linux64DN8FDFolderProfile.pubxml
mkdir -p Linux64DN8FDU20
cp BuildFilter/bin/Release/net8.0/publish/linux-x64/BuildFilter ./Linux64DN8FDU20/
