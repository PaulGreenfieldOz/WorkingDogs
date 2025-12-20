dotnet publish ./Clancy.sln -c release /p:PublishProfile=Linux64DN8FDFolderProfile.pubxml
mkdir -p Linux64DN8FDU24
cp Clancy/bin/Release/net8.0/publish/linux-x64/Clancy ./Linux64DN8FDU24/
