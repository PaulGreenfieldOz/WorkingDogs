dotnet publish .\Pup.sln -c release /p:PublishProfile=Win64DN10AOTFolderProfile.pubxml
mkdir -force Win64DN10AOT
copy Pup\bin\Release\net10.0\publish\win-x64\Pup.exe .\Win64DN10AOT\
