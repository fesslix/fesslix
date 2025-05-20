Write-Output 'Cloning vcpkg repository...'
git clone https://github.com/Microsoft/vcpkg.git
.\vcpkg\bootstrap-vcpkg.bat
.\vcpkg\vcpkg integrate install

Write-Output 'Installing dependencies ...'
.\vcpkg\vcpkg install gsl:x64-windows-static-md gsl:x86-windows-static-md gsl:arm64-windows-static-md boost-format:x64-windows-static-md boost-math:x64-windows-static-md boost-concept-check:x64-windows-static-md boost-random:x64-windows-static-md boost-algorithm:x64-windows-static-md nlopt:x64-windows-static-md nlopt:x86-windows-static-md nlopt:arm64-windows-static-md

.\vcpkg\vcpkg integrate install

Write-Output 'Check if pkg-config Can Find GSL ...'
.\vcpkg\vcpkg list | Select-String gsl
.\vcpkg\vcpkg list | Select-String boost
$env:PKG_CONFIG_PATH="D:\a\fesslix\fesslix\vcpkg\installed\x64-windows-static-md\lib\pkgconfig"
Write-Output '   pkg-config ...'
pkg-config --cflags gsl
pkg-config --libs gsl
pkg-config --modversion gsl

Write-Output 'List files in directory'
#Get-ChildItem -Path "D:\a\fesslix\fesslix\vcpkg" -Recurse -File | Select-Object -ExpandProperty FullName
#dumpbin /symbols D:\a\fesslix\fesslix\vcpkg\installed\x86-windows-static\lib\gsl.lib
#llvm-objdump -p D:\a\fesslix\fesslix\vcpkg\installed\x86-windows-static\lib\gsl.lib
#llvm-readobj D:\a\fesslix\fesslix\vcpkg\installed\x86-windows-static\lib\gsl.lib

