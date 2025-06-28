dotnet publish ./ExtractHeadersFromFasta.csproj -c release /p:PublishProfile=Linux64DN8AOTFolderProfile.pubxml
mkdir -p Linux64DN8AOTU20
cp bin/Release/net8.0/publish/linux-x64/ExtractHeadersFromFasta ./Linux64DN8AOTU20/
