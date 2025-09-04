dotnet publish .\Pup.sln -c release /p:PublishProfile=Win64DN8AOTFolderProfile.pubxml
mkdir -force Win64DN8AOT
copy Pup\bin\Release\net8.0\publish\win-x64\Pup.exe .\Win64DN8AOT\
