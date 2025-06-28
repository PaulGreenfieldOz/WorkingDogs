dotnet publish ./CondenseProkkaTbl.csproj -c release /p:PublishProfile=Linux64DN7FDFolderProfile.pubxml
mkdir -p Linux64DN7FD
cp bin/Release/net7.0/publish/linux-x64/CondenseProkkaTbl ./Linux64DN7FD/
