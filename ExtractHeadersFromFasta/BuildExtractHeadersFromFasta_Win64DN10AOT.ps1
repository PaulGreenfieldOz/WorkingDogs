dotnet publish .\ExtractHeadersFromFasta.sln -c release /p:PublishProfile=Win64DN10AOTFolderProfile.pubxml
mkdir -force Win64DN10AOT
copy ExtractHeadersFromFasta\bin\Release\net10.0\publish\win-x64\ExtractHeadersFromFasta.exe .\Win64DN10AOT\
