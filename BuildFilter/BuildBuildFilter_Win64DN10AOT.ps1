dotnet publish .\BuildFilter.sln -c release /p:PublishProfile=Win64DN10AOTFolderProfile.pubxml
mkdir -force Win64DN10AOT
copy BuildFilter\bin\Release\net10.0\publish\win-x64\BuildFilter.exe .\Win64DN10AOT\
