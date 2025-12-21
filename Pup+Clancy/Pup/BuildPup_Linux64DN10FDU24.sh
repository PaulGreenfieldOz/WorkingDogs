dotnet publish ./Pup.sln -c release /p:PublishProfile=Linux64DN10FDFolderProfile.pubxml
mkdir -p Linux64DN10FDU24
cp Pup/bin/Release/net10.0/publish/linux-x64/Pup ./Linux64DN10FDU24/
