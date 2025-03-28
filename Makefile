# This Makefile is used to compile the FlashRMSD and DockRMSD programs

all: FlashRMSD DockRMSDExt DockRMSD

DockRMSDExt:
	gcc tools/src/DockRMSD/dock_rmsd_modified.c -o tools/DockRMSDExt -O3 -lm

DockRMSD:
	gcc tools/src/DockRMSD/dock_rmsd_original.c -o tools/DockRMSD -O3 -lm

FlashRMSD:
	gcc -I./include ./src/*.c -o tools/FlashRMSD -O3 -lm
