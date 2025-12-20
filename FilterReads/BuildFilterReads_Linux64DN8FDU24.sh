dotnet publish ./FilterReads.sln -c release /p:PublishProfile=Linux64DN8FDFolderProfile.pubxml
mkdir -p Linux64DN8FDU24
cp FilterReads/bin/Release/net8.0/publish/linux-x64/FilterReads ./Linux64DN8FDU24/
