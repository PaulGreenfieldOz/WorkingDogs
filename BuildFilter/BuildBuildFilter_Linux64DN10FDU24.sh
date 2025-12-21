dotnet publish ./BuildFilter.sln -c release /p:PublishProfile=Linux64DN10FDFolderProfile.pubxml
mkdir -p Linux64DN10FDU24
cp BuildFilter/bin/Release/net10.0/publish/linux-x64/BuildFilter ./Linux64DN10FDU24/
