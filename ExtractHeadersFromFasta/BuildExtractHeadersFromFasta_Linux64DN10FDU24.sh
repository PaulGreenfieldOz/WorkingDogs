dotnet publish ./ExtractHeadersFromFasta.sln -c release /p:PublishProfile=Linux64DN10FDFolderProfile.pubxml
mkdir -p Linux64DN10FDU24
cp ExtractHeadersFromFasta/bin/Release/net10.0/publish/linux-x64/ExtractHeadersFromFasta ./Linux64DN10FDU24/
