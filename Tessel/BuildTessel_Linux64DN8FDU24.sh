dotnet publish ./Tessel.sln -c release /p:PublishProfile=Linux64DN8FDFolderProfile.pubxml
mkdir -p Linux64DN8FDU24
cp Tessel/bin/Release/net8.0/publish/linux-x64/Tessel ./Linux64DN8FDU24/
