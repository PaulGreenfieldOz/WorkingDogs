dotnet publish ./Blue.csproj -c release /p:PublishProfile=Linux64DN7AOTFolderProfile.pubxml
mkdir Linux64DN7AOT
cp bin/Release/net7.0/publish/linux-x64/Blue ./Linux64DN7AOT/
