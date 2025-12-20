dotnet publish .\CondenseProkkaTbl.sln -c release /p:PublishProfile=Win64DN10AOTFolderProfile.pubxml
mkdir -force Win64DN10AOT
copy CondenseProkkaTbl\bin\Release\net10.0\publish\win-x64\CondenseProkkaTbl.exe .\Win64DN10AOT\
