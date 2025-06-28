dotnet publish ./CondenseProkkaTbl.csproj -c release /p:PublishProfile=Linux64DN8FDFolderProfile.pubxml
mkdir -p Linux64DN8FD
cp bin/Release/net8.0/publish/linux-x64/CondenseProkkaTbl ./Linux64DN8FD/
