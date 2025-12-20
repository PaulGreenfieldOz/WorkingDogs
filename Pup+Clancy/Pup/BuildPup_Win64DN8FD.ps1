dotnet publish .\Pup.sln -c release /p:PublishProfile=Win64DN8FDFolderProfile.pubxml
mkdir -force Win64DN8FD
copy Pup\bin\Release\net8.0\publish\win-x64\Pup.exe .\Win64DN8FD\

