dotnet publish .\Kelpie_v2.sln -c release /p:PublishProfile=Win64DN8FDFolderProfile.pubxml
mkdir -force Win64DN8FD
copy Kelpie_v2\bin\Release\net8.0\publish\win-x64\Kelpie_v2.exe .\Win64DN8FD\
copy Kelpie_v2\bin\Release\net8.0\publish\win-x64\Kelpie_v2.exe .\Win64DN8FD\Kelpie.exe
