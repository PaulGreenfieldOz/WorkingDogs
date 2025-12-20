dotnet publish .\Tessel.sln -c release /p:PublishProfile=Win64DN8AOTFolderProfile.pubxml
mkdir -force Win64DN8AOT
copy Tessel\bin\Release\net8.0\publish\win-x64\Tessel.exe .\Win64DN8AOT\
