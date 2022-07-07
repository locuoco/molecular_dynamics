#!/bin/bash

# If using WSL, add /mnt to PRUNEPATHS in /etc/updatedb.conf in order to avoid indexing Windows files.
sudo apt install g++-10 libglew-dev libglfw3-dev libfreetype6-dev
sudo update-alternatives --remove-all g++
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 20