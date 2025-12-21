dotnet publish ./Tessel.sln -c release /p:PublishProfile=Linux64DN10FDFolderProfile.pubxml
mkdir -p Linux64DN10FDU24
cp Tessel/bin/Release/net10.0/publish/linux-x64/Tessel ./Linux64DN10FDU24/
