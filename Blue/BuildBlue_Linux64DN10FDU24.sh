dotnet publish ./Blue.sln -c release /p:PublishProfile=Linux64DN10FDFolderProfile.pubxml
mkdir -p Linux64DN10FDU24
cp Blue/bin/Release/net10.0/publish/linux-x64/Blue ./Linux64DN10FDU24/
