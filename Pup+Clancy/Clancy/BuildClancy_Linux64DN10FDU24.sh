dotnet publish ./Clancy.sln -c release /p:PublishProfile=Linux64DN10FDFolderProfile.pubxml
mkdir -p Linux64DN10FDU24
cp Clancy/bin/Release/net10.0/publish/linux-x64/Clancy ./Linux64DN10FDU24/
