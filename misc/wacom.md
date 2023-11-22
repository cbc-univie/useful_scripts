Installing Wacom on Linux
=========================

Go to https://linuxwacom.github.io/ and select Kernel Driver, which brings you to https://github.com/linuxwacom/input-wacom/wiki/Installing-input-wacom-from-source. Install Prerequisities as stated, then proceed with Download (https://github.com/linuxwacom/input-wacom/releases). I used the "Source code (zip)" option.
unzip
go into folder (save somewhere in scratch to be able to access as root)
./autogen.sh
make
as root user: sudo make install
reboot


Problem: only mirror option is available by default
enable as second(third) screen with (https://askubuntu.com/questions/839161/limit-a-graphics-tablet-to-one-monitor):
xinput map-to-output 8 HDMI-0
where 8 is the stylus id (get with xinput (install with sudo apt-get install xinput))
get display with:
xrandr -> in my case HMDI-0 did not work but "HEAD-1" as stated in the "NVIDIA binary driver" section (https://github.com/linuxwacom/xf86-input-wacom/wiki/Dual-and-Multi-Monitor-Set-Up)


Another option, but this one worked only with two but not three screens in my case
(https://github.com/linuxwacom/xf86-input-wacom/wiki/Dual-and-Multi-Monitor-Set-Up):
xsetwacom set "Wacom One Pen Display 13 Pen stylus" MapToOutput "HEAD-1"
get name of tablet with:
xinput; in our case: Wacom One Pen Display 13 Pen stylus
get display with:
xrandr -> using this command HMDI-0 did not work but "HEAD-1" as stated in the "NVIDIA binary driver" section (https://github.com/linuxwacom/xf86-input-wacom/wiki/Dual-and-Multi-Monitor-Set-Up)


Screen Sharing with BigBlueButton:
Screen sharing within BBB for an application opened on the tablet did not work, just black screen appeared.
Solution (having 2 screens, one of them is mirrored with the WACOM): open "NVIDIA X Server Settings", go to "X Server Display Configuration" and select a monitor and for its position, select "Same as" and then the name of your monitor. (https://askubuntu.com/questions/576870/how-can-i-mirror-one-of-the-screens-in-a-3-monitor-setup)
