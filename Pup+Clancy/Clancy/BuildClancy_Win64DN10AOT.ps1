dotnet publish .\Clancy.sln -c release /p:PublishProfile=Win64DN10AOTFolderProfile.pubxml
mkdir -force Win64DN10AOT
copy Clancy\bin\Release\net10.0\publish\win-x64\Clancy.exe .\Win64DN10AOT\
