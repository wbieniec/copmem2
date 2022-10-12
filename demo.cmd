@echo off
REM Demo for copmem2 for Windows
REM Please compile under Visual Studio as Release
REM Plase provide valid dataset files

set R=..\datasets\hum.all.fa
set Q=..\datasets\panTro3.fa
set tst=hp
for %%i in (200 100 80) do (
  for %%j in (8 4 1) do (
    for %%k in (a b c) do (
      x64\Release\MemoryFill 28G
      x64\Release\copmem2.exe  -t %%j -l %%i -o %tst%-%%i-t%%j.txt %R% %Q%
      del %tst%-%%i-t%%j.txt

      x64\Release\MemoryFill 28G
      x64\Release\copmem2.exe  -mf -t %%j -l %%i -o %tst%-%%i-t%%j.txt %R% %Q%
      del %tst%-%%i-t%%j.txt
    )
  )
)