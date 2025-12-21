dotnet publish ./ExtractHeadersFromFasta.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml
mkdir -p Linux64DN10AOTU24
cp ExtractHeadersFromFasta/bin/Release/net10.0/publish/linux-x64/ExtractHeadersFromFasta ./Linux64DN10AOTU24/
