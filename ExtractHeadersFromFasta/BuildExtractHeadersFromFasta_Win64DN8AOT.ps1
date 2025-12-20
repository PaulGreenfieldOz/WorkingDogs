dotnet publish .\ExtractHeadersFromFasta.sln -c release /p:PublishProfile=Win64DN8AOTFolderProfile.pubxml
mkdir -force Win64DN8AOT
copy ExtractHeadersFromFasta\bin\Release\net8.0\publish\win-x64\ExtractHeadersFromFasta.exe .\Win64DN8AOT\
