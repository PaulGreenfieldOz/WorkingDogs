dotnet publish .\CondenseProkkaTbl.csproj -c release /p:PublishProfile=Win64DN8AOTFolderProfile.pubxml
mkdir -force Win64DN8AOT
copy bin\Release\net8.0\publish\win-x64\CondenseProkkaTbl.exe .\Win64DN8AOT\
