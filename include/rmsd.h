#ifndef RMSD_H
#define RMSD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "common.h"

#include "utils.h"
#include "io.h"

int check_for_basic_matches(MolecularData*, MolecularData*);
double rmsd(MolecularData*, MolecularData*, int*);

#endif // RMSD_H