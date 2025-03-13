# This Makefile is used to compile the Frmsd and DockRMSD programs

all: Frmsd DockRMSDExt DockRMSD

DockRMSDExt:
	gcc tools/src/DockRMSD/dock_rmsd_modified.c -o tools/DockRMSDExt -O3 -lm

DockRMSD:
	gcc tools/src/DockRMSD/dock_rmsd_original.c -o tools/DockRMSD -O3 -lm

Frmsd:
	gcc -I./include ./src/*.c -o tools/Frmsd -O3 -lm
