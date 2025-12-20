dotnet publish .\Tessel.sln -c release /p:PublishProfile=Win64DN10FDFolderProfile.pubxml
mkdir -force Win64DN10FD
copy Tessel\bin\Release\net10.0\publish\win-x64\Tessel.exe .\Win64DN10FD\
