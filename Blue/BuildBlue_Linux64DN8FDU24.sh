dotnet publish ./Blue.sln -c release /p:PublishProfile=Linux64DN8FDFolderProfile.pubxml
mkdir -p Linux64DN8FDU24
cp Blue/bin/Release/net8.0/publish/linux-x64/Blue ./Linux64DN8FDU24/
