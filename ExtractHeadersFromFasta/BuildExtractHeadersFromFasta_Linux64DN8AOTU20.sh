dotnet publish ./ExtractHeadersFromFasta.sln -c release /p:PublishProfile=Linux64DN8AOTFolderProfile.pubxml
mkdir -p Linux64DN8AOTU20
cp ExtractHeadersFromFasta/bin/Release/net8.0/publish/linux-x64/ExtractHeadersFromFasta ./Linux64DN8AOTU20/
