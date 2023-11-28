dotnet publish .\Kelpie_v2.csproj -c release /p:PublishProfile=Win64DN8AOTFolderProfile.pubxml
mkdir -force Win64DN8AOT
copy bin\Release\net8.0\publish\win-x64\Kelpie_v2.exe .\Win64DN8AOT\
