dotnet publish ./CondenseProkkaTbl.sln -c release /p:PublishProfile=Linux64DN10FDFolderProfile.pubxml
mkdir -p Linux64DN10FDU24
cp CondenseProkkaTbl/bin/Release/net10.0/publish/linux-x64/CondenseProkkaTbl ./Linux64DN10FDU24/
