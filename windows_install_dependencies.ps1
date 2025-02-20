Write-Output 'Cloning vcpkg repository...'
git clone https://github.com/microsoft/vcpkg.git
Write-Output 'Bootstrapping vcpkg...'
.\vcpkg\bootstrap-vcpkg.bat
Write-Output 'Installing GSL and Boost...'
.\vcpkg\vcpkg install gsl boost
Write-Output 'Setting environment variables...'
[System.Environment]::SetEnvironmentVariable('VCPKG_ROOT', "$PWD\vcpkg", [System.EnvironmentVariableTarget]::Process)
[System.Environment]::SetEnvironmentVariable('CMAKE_TOOLCHAIN_FILE', "$PWD\vcpkg\scripts\buildsystems\vcpkg.cmake", [System.EnvironmentVariableTarget]::Process)
Write-Output 'Setup completed successfully.'

