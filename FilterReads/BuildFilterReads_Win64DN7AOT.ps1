dotnet publish .\FilterReads.csproj -c release /p:PublishProfile=Win64DN7AOTFolderProfile.pubxml
mkdir -force Win64DN7AOT
copy bin\Release\net7.0\publish\win-x64\FilterReads.exe .\Win64DN7AOT\
