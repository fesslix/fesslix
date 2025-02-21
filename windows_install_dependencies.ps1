Write-Output 'Cloning vcpkg repository...'
#git clone https://github.com/Microsoft/vcpkg.git
#.\vcpkg\bootstrap-vcpkg.bat
vcpkg integrate install

Write-Output 'Installing dependencies ...'
vcpkg install gsl gsl:x86-windows boost-format boost-math boost-concept-check boost-random

vcpkg integrate install

#git clone https://github.com/microsoft/vcpkg.git
#Write-Output 'Bootstrapping vcpkg...'
#.\vcpkg\bootstrap-vcpkg.bat
#Write-Output 'Installing GSL and Boost...'
#.\vcpkg\vcpkg install gsl boost
#Write-Output 'Setup completed successfully.'

#Write-Output 'VCPKG update an version'
#vcpkg version
#vcpkg update
#vcpkg upgrade --no-dry-run
#vcpkg version

#Write-Output 'Installing dependencies ...'
#vcpkg install --triplet x64-windows

Write-Output 'Check if pkg-config Can Find GSL ...'
vcpkg list | Select-String gsl
$env:PKG_CONFIG_PATH="C:\vcpkg\installed\x64-windows\lib\pkgconfig"
Write-Output '   pkg-config ...'
pkg-config --cflags gsl
pkg-config --libs gsl
pkg-config --modversion gsl

#Write-Output 'Setting environment variables...'
#[System.Environment]::SetEnvironmentVariable('VCPKG_ROOT', "$PWD\vcpkg", [System.EnvironmentVariableTarget]::Process)
#[System.Environment]::SetEnvironmentVariable('CMAKE_TOOLCHAIN_FILE', "$PWD\vcpkg\scripts\buildsystems\vcpkg.cmake", [System.EnvironmentVariableTarget]::Process)
#Write-Output "20250221: $CMAKE_TOOLCHAIN_FILE"

Write-Output 'List files in directory'
#Get-ChildItem -Path "C:\vcpkg" -Recurse -File | Select-Object -ExpandProperty FullName
