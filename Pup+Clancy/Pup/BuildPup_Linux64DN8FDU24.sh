dotnet publish ./Pup.sln -c release /p:PublishProfile=Linux64DN8FDFolderProfile.pubxml
mkdir -p Linux64DN8FDU24
cp Pup/bin/Release/net8.0/publish/linux-x64/Pup ./Linux64DN8FDU24/
