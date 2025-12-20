dotnet publish .\Clancy.sln -c release /p:PublishProfile=Win64DN10FDFolderProfile.pubxml
mkdir -force Win64DN10FD
copy Clancy\bin\Release\net10.0\publish\win-x64\Clancy.exe .\Win64DN10FD\

