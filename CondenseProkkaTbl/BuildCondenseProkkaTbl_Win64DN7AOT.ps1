dotnet publish .\CondenseProkkaTbl.csproj -c release /p:PublishProfile=Win64DN7AOTFolderProfile.pubxml
mkdir -force Win64DN7AOT
copy bin\Release\net7.0\publish\win-x64\CondenseProkkaTbl.exe .\Win64DN7AOT\
