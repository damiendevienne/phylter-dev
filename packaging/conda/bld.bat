@echo off
cmake -S . -B build-conda -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" -DPHYLTER2_OPENMP=OFF
if errorlevel 1 exit /b 1
cmake --build build-conda --parallel %CPU_COUNT%
if errorlevel 1 exit /b 1
ctest --test-dir build-conda --output-on-failure
if errorlevel 1 exit /b 1
cmake --install build-conda --component Runtime
if errorlevel 1 exit /b 1
