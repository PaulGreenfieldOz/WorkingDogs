dotnet publish .\FilterReads.sln -c release /p:PublishProfile=Win64DN10AOTFolderProfile.pubxml
mkdir -force Win64DN10AOT
copy FilterReads\bin\Release\net10.0\publish\win-x64\FilterReads.exe .\Win64DN10AOT\
