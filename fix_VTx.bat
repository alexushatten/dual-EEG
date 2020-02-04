:: This Windows batch script is useful on some PCs with Windows 10 for which VirtualBox does not function properly (anymore)

:: NB: I had an error saying that VTx was not enabled when starting VirtualBox, and despite trying all possible fixes
:: mentioned on VirtualBox community forum (including ticking and unticking many boxes (Hyper-V, etc) in Windows features),
:: I realized it (also) had to do with DeviceGuard. DG Readiness Tool helped me to fix it when issuing the following commands:

:: Make sure command prompt's current directory is in proper drive
C:
:: Disable Virtualization Based Security
reg add HKLM\System\CurrentControlSet\Control\DeviceGuard /t REG_DWORD /v EnableVirtualizationBasedSecurity /d 0 /f
:: Specify the platform security level at the next reboot: turn on VBS with Secure Boot and direct memory access (DMA)
reg add HKLM\System\CurrentControlSet\Control\DeviceGuard /t REG_DWORD /v RequirePlatformSecurityFeatures /d 3 /f
:: Disable Credential Guard
reg add HKLM\System\CurrentControlSet\Control\Lsa /t REG_DWORD /v LsaCfgFlags /d 0 /f
:: Run Device Guard and Credential Guard hardware readiness tool, disable Device Guard and Credential Guard, and restart PC
powershell -executionpolicy bypass ""C:\Temp\dgreadiness_v3.6\DG_Readiness_Tool_v3.6.ps1 -Disable -AutoReboot""

:: This script will restart your computer. After every Windows update, I had to run it again, likely due to institutional policies.
