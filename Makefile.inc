# PYTHIA configuration file.
# Generated on Fri 24 Jan 11:20:36 -05 2020 with the user supplied options:
# --enable-shared
# --with-root=/home/sane/builds/root/

# Install directory prefixes.
PREFIX_BIN=/home/sane/builds/pythia8244/bin
PREFIX_INCLUDE=/home/sane/builds/pythia8244/include
PREFIX_LIB=/home/sane/builds/pythia8244/lib
PREFIX_SHARE=/home/sane/builds/pythia8244/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=true
CXX=g++
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

# ROOT configuration.
ROOT_USE = true
ROOT_BIN = /home/sane/builds/root/bin/
ROOT_INCLUDE = /home/sane/builds/root/include/
ROOT_LIB = /home/sane/builds/root/lib/