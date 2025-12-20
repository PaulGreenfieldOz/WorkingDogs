dotnet publish .\Kelpie_v2.sln -c release /p:PublishProfile=Win64DN10AOTFolderProfile.pubxml
mkdir -force Win64DN10AOT
copy Kelpie_v2\bin\Release\net10.0\publish\win-x64\Kelpie_v2.exe .\Win64DN10AOT\
copy Kelpie_v2\bin\Release\net10.0\publish\win-x64\Kelpie_v2.exe .\Win64DN10AOT\Kelpie.exe
