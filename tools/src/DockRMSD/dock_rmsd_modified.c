#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/*
 DockRMSD: an open-source tool for atom mapping and RMSD calculation of symmetric molecules through graph isomorphism

 Written by Eric Bell
 v1.0 written 5/2/2019
 latest update (v1.1) written 8/26/2019

 To compile, use:
 gcc DockRMSD.c -o DockRMSD -lm -O3
*/

const int MAXBONDS = 15;        // Maximum number of bonds allowable on a single atom
const int MAXLINELENGTH = 150; // Maximum length (in characters) of a line in a mol2 file
const int MAXBONDSTRING = 32;   // Maximum length (in characters) of a bond string in a mol2 file
const double MAXMAPCOUNT = 0;  // Maximum amount of possible mappings before symmetry heuristic is used
const int MAXDEPTH = 2;
const char *header = "######################################################################\n"
                     "# DockRMSD (v1.1): docking pose distance calculation                 #\n"
                     "#                                                                    #\n"
                     "# If you use DockRMSD in your work, please cite:                     #\n"
                     "#                                                                    #\n"
                     "# Bell, E.W., Zhang, Y. DockRMSD: an open-source tool for atom       #\n"
                     "# mapping and RMSD calculation of symmetric molecules through graph  #\n"
                     "# isomorphism. Journal of Cheminformatics, 11:40 (2019).             #\n"
                     "######################################################################\n";

int grabAtomCount(FILE *mol2, int hflag);
int arrayIdentity(char **arr1, char **arr2, int arrlen);
int inArray(int n, int *arr, int arrlen);
int readMol2Blocks(char ****atoms, double ****coords, char *****bonds, int ***nums, int **atom_counts, int **bond_counts, FILE *mol2, int hflag);
void readMol2(char **atoms, double **coords, char ***bonds, int *nums, FILE *mol2, int atomcount, int hflag);
int generalizeBonds(char ***bonds, int atomcount);
char **buildTree(int depth, int index, char **atoms, char ***bonds, char *prestring, int prevind, int atomcount);
double searchAssigns(int atomcount, int **allcands, int candcounts[], int *assign, char ***tempbond, char ***querybond, double **querycoord, double **tempcoord, int *bestassign);
double assignAtoms(char **tempatom, char ***tempbond, char **queryatom, char ***querybond, double **querycoord, double **tempcoord, int *querynums, int *tempnums, int atomcount, int simpleflag);
int validateBonds(int *atomassign, int proposedatom, int assignpos, char ***querybond, char ***tempbond, double querydists[], int atomcount);

int main(int argc, char *argv[])
{
  // Check if all inputs are given and valid
  if (argc < 3 || argc > 7)
  {
    printf("Program must have a query and template mol2 filename, and at most 3 options\n");
    return 1;
  }

  int flags_start_id = 2;
  FILE *template = fopen(argv[1], "r");
  FILE *query = NULL;
  if (argv[2][0] != '-')
  {
    query = fopen(argv[2], "r");
    flags_start_id = 3;
  }

  // Runtime flag processing
  int hflag = 0;
  int corrflag = 0;
  int simpleflag = 0;
  int crossflag = 0;
  int i;
  for (i = flags_start_id; i < argc; i++)
  {
    if (!strcmp(argv[i], "-h"))
    {
      hflag = 1;
    }
    else if (!strcmp(argv[i], "-c"))
    {
      corrflag = 1;
    }
    else if (!strcmp(argv[i], "-s"))
    {
      simpleflag = 1;
    }
    else if (!strcmp(argv[i], "-x"))
    {
      crossflag = 1;
    }
    else
    {
      printf("Unrecognized flag %s. Exiting...\n", argv[i]);
      return 5;
    }
  }

  if (template == NULL)
  {
    printf("Template file %s not found!\n", argv[2]);
    return 1;
  }
  if (query == NULL && !crossflag)
  {
    printf("Query file %s not found!\n", argv[1]);
    return 1;
  }

  if (!simpleflag)
  {
    printf("%s\n", header);
  }

  if (crossflag) {
    char ***atoms;
    double ***coords;
    char ****bonds;
    int *atom_counts;
    int *bond_counts;
    int **nums;

    int blocks_count = readMol2Blocks(&atoms, &coords, &bonds, &nums, &atom_counts, &bond_counts, template, hflag);

    for (int i = 0; i < blocks_count; i++) {
      generalizeBonds(bonds[i], atom_counts[i]);
    }


    double **best_rmsds = malloc(sizeof(double *) * blocks_count);
    for (int i = 0; i < blocks_count; i++) {
      best_rmsds[i] = malloc(sizeof(double) * blocks_count);
      for (int j = 0; j < blocks_count; j++) {
        best_rmsds[i][j] = DBL_MAX;
      }

      for (int j = 0; j < blocks_count; j++) {
        if (i > j) {
          best_rmsds[i][j] = best_rmsds[j][i];
        } else if (i == j) {
          best_rmsds[i][j] = 0.0;
        } else {
          char **queryatoms = atoms[i];
          double **querycoords = coords[i];
          char ***querybonds = bonds[i];
          int *querynums = nums[i];
          int querycount = atom_counts[i];
          
          char **tempatoms = atoms[j];
          double **tempcoords = coords[j];
          char ***tempbonds = bonds[j];
          int *tempnums = nums[j];
          int tempcount = atom_counts[j];

          char **flatquerybonds = malloc(sizeof(char *) * (querycount * querycount));
          char **flattempbonds = malloc(sizeof(char *) * (tempcount * tempcount));
          
          for (int k = 0; k < querycount; k++)
          {
            memcpy(flatquerybonds + querycount * k, *(querybonds + k), sizeof(char *) * querycount);
            memcpy(flattempbonds + tempcount * k, *(tempbonds + k), sizeof(char *) * tempcount);
          }

          // Find all mappings of template atoms onto query atoms and count how many valid mappings exist
          best_rmsds[i][j] = assignAtoms(tempatoms, tempbonds, queryatoms, querybonds, querycoords, tempcoords, querynums, tempnums, querycount, simpleflag);
        }
      }
    }

    for (int i = 0; i < blocks_count; i++) {
      for (int j = 0; j < blocks_count; j++) {
        printf("%.3f ", best_rmsds[i][j]);
      }
      printf("\n");
    }
  }
  else {
    // Count how many atoms are in the query and template
    int querycount = grabAtomCount(query, hflag);
    int tempcount = grabAtomCount(template, hflag);
    if (querycount != tempcount)
    {
      printf("%s\n", "Error: Query and template don't have the same atom count!");
      return 2;
    }

    if (querycount == 0)
    {
      printf("%s\n", "Error: Query file has no atoms!");
      return 2;
    }
    if (tempcount == 0)
    {
      printf("%s\n", "Error: Template file has no atoms!");
      return 2;
    }

    // Initialize pointer arrays and fill them with information from the mol2 files
    char **queryatoms = malloc(sizeof(char *) * querycount);
    double **querycoords = malloc(sizeof(double *) * querycount);
    char ***querybonds = malloc(sizeof(char **) * querycount);
    char **tempatoms = malloc(sizeof(char *) * tempcount);
    double **tempcoords = malloc(sizeof(double *) * tempcount);
    char ***tempbonds = malloc(sizeof(char **) * tempcount);
    int *querynums = malloc(sizeof(int) * querycount);
    int *tempnums = malloc(sizeof(int) * tempcount);
    int j;
    for (i = 0; i < querycount; i++)
    {
      char *queryatom = malloc(sizeof(char) * 3);
      *(queryatoms + i) = queryatom;
      char *tempatom = malloc(sizeof(char) * 3);
      *(tempatoms + i) = tempatom;
      double *querycoord = malloc(sizeof(double) * 3);
      *(querycoords + i) = querycoord;
      double *tempcoord = malloc(sizeof(double) * 3);
      *(tempcoords + i) = tempcoord;
      char **querybondrow = malloc(sizeof(char *) * querycount);
      char **tempbondrow = malloc(sizeof(char *) * tempcount);
      for (j = 0; j < querycount; j++)
      {
        char *querybond = malloc(sizeof(char) * 3);
        strcpy(querybond, "");
        *(querybondrow + j) = querybond;
        char *tempbond = malloc(sizeof(char) * 3);
        strcpy(tempbond, "");
        *(tempbondrow + j) = tempbond;
      }
      *(querybonds + i) = querybondrow;
      *(tempbonds + i) = tempbondrow;
    }
    readMol2(queryatoms, querycoords, querybonds, querynums, query, querycount, hflag);
    readMol2(tempatoms, tempcoords, tempbonds, tempnums, template, tempcount, hflag);
    fclose(query);
    fclose(template);
    generalizeBonds(querybonds,querycount);
    generalizeBonds(tempbonds,tempcount);
    // Check that the query and template are in fact the same molecules
    if (!arrayIdentity(queryatoms, tempatoms, querycount))
    {
      printf("%s\n", "Template and query don't have the same atoms.");
      return 3;
    }

    char **flatquerybonds = malloc(sizeof(char *) * (querycount * querycount));
    char **flattempbonds = malloc(sizeof(char *) * (tempcount * tempcount));
    for (i = 0; i < querycount; i++)
    {
      memcpy(flatquerybonds + querycount * i, *(querybonds + i), sizeof(char *) * querycount);
      memcpy(flattempbonds + tempcount * i, *(tempbonds + i), sizeof(char *) * tempcount);
    }
    if (!arrayIdentity(flatquerybonds, flattempbonds, querycount * querycount))
    {
      // Remove bond typing if they don't agree between query and template
      generalizeBonds(querybonds, querycount);
      generalizeBonds(tempbonds, tempcount);
      for (i = 0; i < querycount; i++)
      {
        memcpy(flatquerybonds + querycount * i, *(querybonds + i), sizeof(char *) * querycount);
        memcpy(flattempbonds + tempcount * i, *(tempbonds + i), sizeof(char *) * tempcount);
      }
      // If the general bonds still don't agree, the molecules aren't the same
      if (!arrayIdentity(flatquerybonds, flattempbonds, querycount * querycount))
      {
        printf("%s\n", "Template and query don't have the same bonding network.");
        return 3;
      }
    }
    free(flatquerybonds);
    free(flattempbonds);

    // Find all mappings of template atoms onto query atoms and count how many valid mappings exist
    if (!corrflag)
    {
      printf("%.3f\n", assignAtoms(tempatoms, tempbonds, queryatoms, querybonds, querycoords, tempcoords, querynums, tempnums, querycount, simpleflag));
    }
    else
    {
      for (i = 0; i < querycount; i++)
      {
        if (strcmp(*(queryatoms + i), *(tempatoms + i)))
        {
          printf("Atomic correspondence assumption is not valid: atom %d is %s in query, %s in template.\n", i, *(queryatoms + i), *(tempatoms + i));
          return 5;
        }
      }
      double squaredev = 0.0;
      for (i = 0; i < querycount; i++)
      {
        for (j = 0; j < 3; j++)
        {
          squaredev += pow(*(*(querycoords + i) + j) - *(*(tempcoords + i) + j), 2.0);
        }
      }
      double dockrmsd = pow(squaredev / ((double)querycount), 0.5);
      if (simpleflag)
      {
        printf("%.3f\n", dockrmsd);
      }
      else
      {
        printf("Calculated Docking RMSD: %.3f\n", dockrmsd);
      }
      }
    }
  return 0;
}

// Comparator for compatibility with qsort
int strcompar(const void *a, const void *b) { return strcmp(*(char **)a, *(char **)b); }

// Returns 1 if two arrays contain the same string elements, otherwise returns 0
int arrayIdentity(char **arr1, char **arr2, int arrlen)
{
  char *list1[arrlen];
  char *list2[arrlen];
  int i;
  for (i = 0; i < arrlen; i++)
  {
    list1[i] = *(arr1 + i);
    list2[i] = *(arr2 + i);
  }
  qsort(list1, arrlen, sizeof(list1[0]), strcompar);
  qsort(list2, arrlen, sizeof(list2[0]), strcompar);
  for (i = 0; i < arrlen; i++)
  {
    if (strcmp(list1[i], list2[i]))
    {
      return 0;
    }
  }
  return 1;
}

// Returns the index+1 if the element n is already in the array, otherwise returns 0
int inArray(int n, int *arr, int arrlen)
{
  int i;
  for (i = 0; i < arrlen; i++)
  {
    if (*(arr + i) == n)
    {
      return i + 1;
    }
  }
  return 0;
}

// Returns the count of atoms in a mol2 file
int grabAtomCount(FILE *mol2, int hflag)
{
  char line[MAXLINELENGTH];
  int atomcount = 0;
  int countflag = 0;
  while (fgets(line, MAXLINELENGTH, mol2) != NULL)
  {
    if (strlen(line) > 1 && line[strlen(line) - 2] == '\r')
    { // For windows line endings
      line[strlen(line) - 2] = '\n';
      line[strlen(line) - 1] = '\0';
    }
    if (!strcmp(line, "@<TRIPOS>ATOM\n"))
    {
      countflag = 1;
      continue;
    }
    if (!strncmp(line, "@<TRIPOS>", 9) && countflag)
    { // If we reach a new section and have already read atoms, break
      countflag = 0;
      break;
    }
    if (countflag && strlen(line) > 1)
    {
      char *token = strtok(line, " \t");
      int i;
      for (i = 0; i < 5; i++)
      {
        token = strtok(NULL, " \t");
      }
      if (hflag || strcmp(token, "H"))
      {
        atomcount++;
      }
    }
  }
  if (ferror(mol2))
  {
    printf("Error %d while reading in file.\n", ferror(mol2));
  }
  rewind(mol2); // resets the file pointer for use in other functions
  return atomcount;
}

int readMol2Blocks(char ****atoms, double ****coords, char *****bonds, int ***nums, int **atom_counts, int **bond_counts, FILE *mol2, int hflag) {
  char line[MAXLINELENGTH];
  /* These will hold one block per molecule read */
  char ***atomBlocks = NULL;    /* Each element: a char** for one molecule */
  double ***coordBlocks = NULL;   /* Each element: a double** (each row is 3 doubles) */
  char ****bondBlocks = NULL;    /* Each element: a char*** (2D array of bond strings) */
  int **numBlocks = NULL;        /* Each element: an int* (atom numbers) */
  int *atomCounts = NULL;      /* Each element: an int (number of atoms in the block) */
  int *bondCounts = NULL;      /* Each element: an int (number of bonds in the block) */

  int blockCount = 0;
  
  /* We first search for the first block header. */
  while (fgets(line, MAXLINELENGTH, mol2) != NULL) {
      if (strncmp(line, "@<TRIPOS>MOLECULE", 17) == 0)
          break;
  }
  /* If no molecule block is found, set outputs to NULL and return */
  if (feof(mol2)) {
      *atoms = NULL;
      *coords = NULL;
      *bonds = NULL;
      *nums = NULL;
      *atom_counts = NULL;
      *bond_counts = NULL;
      return 0;
  }
  
  /* We use a buffer to “remember” a header line that might belong to the next block. */
  char headerLine[MAXLINELENGTH];
  strcpy(headerLine, line);
  
  /* Process blocks one after another */
  do {
      /* Each new block should begin with "@<TRIPOS>MOLECULE" */
      if (strncmp(headerLine, "@<TRIPOS>MOLECULE", 17) != 0) {
          /* If the line isn’t a molecule header, try to read until one is found */
          if (fgets(headerLine, MAXLINELENGTH, mol2) == NULL)
              break;
          if (strncmp(headerLine, "@<TRIPOS>MOLECULE", 17) != 0)
              continue;
      }
      
      /* Increase block count and reallocate block arrays */
      blockCount++;
      atomBlocks = realloc(atomBlocks, blockCount * sizeof(char **));
      coordBlocks = realloc(coordBlocks, blockCount * sizeof(double **));
      bondBlocks = realloc(bondBlocks, blockCount * sizeof(char ***));
      numBlocks = realloc(numBlocks, blockCount * sizeof(int *));
      atomCounts = realloc(atomCounts, blockCount * sizeof(int));
      bondCounts = realloc(bondCounts, blockCount * sizeof(int));

      /* Read molecule header lines.
       * The first line (molecule name) is read but not used.
       */
      char molName[MAXLINELENGTH];
      if (fgets(molName, MAXLINELENGTH, mol2) == NULL)
          break;
      
      /* The next header line should contain counts.
       * Typically the first number is the atom count.
       */
      char countsLine[MAXLINELENGTH];
      if (fgets(countsLine, MAXLINELENGTH, mol2) == NULL)
          break;
      int atomCount = 0, bondCount = 0;
      sscanf(countsLine, "%d%d", &atomCount, &bondCount);
      
      atomCounts[blockCount - 1] = atomCount;
      bondCounts[blockCount - 1] = bondCount;

      /* Allocate arrays for this block.
       * (For simplicity we allocate atomCount slots even though, if hflag==0,
       *  the number of atoms actually stored may be lower.)
       */
      char **atomsBlock = malloc(atomCount * sizeof(char *));
      double **coordsBlock = malloc(atomCount * sizeof(double *));
      for (int i = 0; i < atomCount; i++)
          coordsBlock[i] = malloc(3 * sizeof(double));
      /* Allocate a 2D array for bonds.
       * We assume a square array of size atomCount x atomCount.
       */
      char ***bondsBlock = malloc(atomCount * sizeof(char **));
      for (int i = 0; i < atomCount; i++) {
          bondsBlock[i] = malloc(atomCount * sizeof(char *));
          for (int j = 0; j < atomCount; j++) {
              bondsBlock[i][j] = calloc(MAXBONDSTRING, sizeof(char));
          }
      }
      int *numsBlock = malloc(atomCount * sizeof(int));

      /* Skip lines until the ATOM section is reached */
      while (fgets(line, MAXLINELENGTH, mol2) != NULL) {
          if (strncmp(line, "@<TRIPOS>ATOM", 13) == 0)
              break;
      }
      
      /* Read the ATOM section.
       * We assume that the atom section continues until a line starting with '@'
       * (which signals either the BOND section or the start of a new block).
       */
      int i = 0;
      int atomNums[atomCount]; /* Temporary array to map original atom numbers */
      while (fgets(line, MAXLINELENGTH, mol2) != NULL) {
          if (line[0] == '@') {
              /* If we see a new section header, check if it is the BOND section. */
              break;
          }
          if (strlen(line) < 2)
              continue; /* Skip blank lines */
          
          /* Example atom line parsing:
           *   field1: atom id (int)
           *   field2: atom name (ignored)
           *   field3-5: x, y, z coordinates (double)
           *   field6: atom type (string) – we use the first letter before any dot.
           */
          double coord[3];
          char *token = strtok(line, " \t");
          int aNum = atoi(token);
          token = strtok(NULL, " \t");  /* Skip atom name */
          for (int j = 0; j < 3; j++) {
              token = strtok(NULL, " \t");
              coord[j] = atof(token);
          }
          token = strtok(NULL, " \t");
          /* If hflag is set or the type is not "H", store the atom. */
          if (hflag || strcmp("H", token) != 0) {
              /* Get element symbol (strip any dot and following characters) */
              char *element = strtok(token, ".");
              atomsBlock[i] = malloc(strlen(element) + 1);
              strcpy(atomsBlock[i], element);
              atomNums[i] = aNum;
              for (int j = 0; j < 3; j++)
                  coordsBlock[i][j] = coord[j];
              i++;
          }
      }

      if (strncmp(line, "@<TRIPOS>BOND", 13) != 0) {
        while (fgets(line, MAXLINELENGTH, mol2) != NULL) {
            if (strncmp(line, "@<TRIPOS>BOND", 13) == 0)
                break;
        }
      }

      /* The actual number of atoms stored */
      int actualAtomCount = i;
      atomCounts[blockCount - 1] = actualAtomCount;
      int file_end = 1;
      /* If the line we last read is the start of the BOND section, process bonds */
      if (strncmp(line, "@<TRIPOS>BOND", 13) == 0) {
          /* Read bond lines until a line starting with '@' is encountered (or EOF) */
          while (fgets(line, MAXLINELENGTH, mol2) != NULL) {
              if (line[0] == '@') {
                  strcpy(headerLine, line);
                  file_end = 0;
                  break;
              }
              /* Parse a bond line.
               * Expected fields: bond id, from atom, to atom, bond type, ...
               */
              char *token = strtok(line, " \t");  /* bond id, ignored */
              token = strtok(NULL, " \t");
              int from = inArray(atoi(token), atomNums, actualAtomCount) - 1;
              token = strtok(NULL, " \t");
              int to = inArray(atoi(token), atomNums, actualAtomCount) - 1;
              token = strtok(NULL, " \t");
              if(token)
                  token = strtok(token, "\n");
              if (from >= 0 && to >= 0) {
                  strcpy(bondsBlock[to][from], token);
                  strcpy(bondsBlock[from][to], token);
              } else {
                bondCounts[blockCount - 1]--;
              }
          }
      } else {
          /* If we did not enter the bond section,
           * then clear headerLine so that we try to read a new block.
           */
          headerLine[0] = '\0';
      }
      
      /* Save atom numbers for this block */
      for (int k = 0; k < actualAtomCount; k++)
          numsBlock[k] = atomNums[k];
      
      /* Store this block’s data in our blocks arrays */
      atomBlocks[blockCount - 1] = atomsBlock;
      coordBlocks[blockCount - 1] = coordsBlock;
      bondBlocks[blockCount - 1] = bondsBlock;
      numBlocks[blockCount - 1] = numsBlock;

      /* If headerLine is empty, try to read next header line from the file */
      if (file_end || headerLine[0] == '\0') {
          if (fgets(headerLine, MAXLINELENGTH, mol2) == NULL)
              break;
      }
  } while (!feof(mol2));

  *atoms = atomBlocks;
  *coords = coordBlocks;
  *bonds = bondBlocks;
  *nums = numBlocks;
  *atom_counts = atomCounts;
  *bond_counts = bondCounts;
  
  /* Store the number of blocks read (accessible globally) */
  return blockCount;
}

// Fills atoms, coords, and bonds with information contained within a mol2 file
void readMol2(char **atoms, double **coords, char ***bonds, int *nums, FILE *mol2, int atomcount, int hflag)
{
  char line[MAXLINELENGTH];
  int atomnums[atomcount]; // Keeps track of all non-H atom numbers for bond reading
  int i = 0;
  int sectionflag = 0; // Value is 1 when reading atoms, 2 when reading bonds, 0 before atoms, >2 after bonds
  int tripos_mol_read = 0;
  while (fgets(line, MAXLINELENGTH, mol2) != NULL)
  {
    if (strlen(line) < 2)
    { // Skip empty lines
      continue;
    }
    if (strlen(line) > 1 && line[strlen(line) - 2] == '\r')
    { // Handling windows line endings
      line[strlen(line) - 2] = '\n';
      line[strlen(line) - 1] = '\0';
    }
    if (!strncmp(line, "@<TRIPOS>MOLECULE", 17))
    { // If we reach a new molecule header, reset the section flag
      if (tripos_mol_read) break;
      tripos_mol_read = 1;
      sectionflag = 0;
      continue;
    }
    else if (!strncmp(line, "@<TRIPOS>ATOM", 13)) 
    {
      sectionflag = 1;
    } 
    else if (!strncmp(line, "@<TRIPOS>BOND", 13))
    {
      sectionflag = 2;
    } 
    else if (!strncmp(line, "@<TRIPOS>", 9))
    { // If we reach a new section and have already read atoms, break
      sectionflag = 0;
    }
    else if (sectionflag == 1)
    { // Reading in atoms and coordinates
      double coord[3];
      int j = 0;
      char *parts = strtok(line, " \t");
      int atomnum = atoi(parts);
      parts = strtok(NULL, " \t");
      for (j = 0; j < 3; j++)
      {
        parts = strtok(NULL, " \t");
        coord[j] = atof(parts);
      }
      parts = strtok(NULL, " \t");
      if (hflag || strcmp("H", parts))
      {
        char *element = strtok(parts, ".");
        strcpy(*(atoms + i), element);
        atomnums[i] = atomnum;
        for (j = 0; j < 3; j++)
        {
          *(*(coords + i) + j) = coord[j];
        }
        i++;
      }
    }
    else if (sectionflag == 2)
    { // Reading in bonding network
      char *parts = strtok(line, " \t");
      parts = strtok(NULL, " \t");
      int from = inArray(atoi(parts), atomnums, atomcount) - 1;
      parts = strtok(NULL, " \t");
      int to = inArray(atoi(parts), atomnums, atomcount) - 1;
      parts = strtok(NULL, " \t");
      parts = strtok(parts, "\n");
      if (from >= 0 && to >= 0)
      {
        strcpy(*(*(bonds + to) + from), parts);
        strcpy(*(*(bonds + from) + to), parts);
      }
    }
  }
  for (i = 0; i < atomcount; i++)
  {
    *(nums + i) = atomnums[i];
  }
}

// Changes all bond types to generic "b" if the bond types don't agree between query and template. Returns true if this has already been done, false if not.
int generalizeBonds(char ***bonds, int atomcount)
{
  int i;
  int j;
  for (i = 0; i < atomcount; i++)
  {
    for (j = 0; j < atomcount; j++)
    {
      if (strcmp(*(*(bonds + i) + j), ""))
      {
        if (!strcmp(*(*(bonds + i) + j), "b"))
        {
          return 1;
        }
        strcpy(*(*(bonds + i) + j), "b");
      }
    }
  }
  return 0;
}

// Recursive function that returns a pointer array of leaves in the bonding tree at a specified depth
char **buildTree(int depth, int index, char **atoms, char ***bonds, char *prestring, int prevind, int atomcount)
{
  int n;
  if (depth == 0)
  { // Base case, if max depth is reached return the prestring
    char **outp = malloc(sizeof(char *) * 2);
    char *outstring = malloc(sizeof(char) * (strlen(prestring) + 1));
    strcpy(outstring, prestring);
    *outp = outstring;
    *(outp + 1) = NULL;
    return outp;
  }
  else
  {
    // Grab all immediate neighbors of the current atom
    int bondinds[MAXBONDS];
    char *bondtypes[MAXBONDS];
    int i;
    int bondi = 0;
    for (i = 0; i < atomcount; i++)
    {
      if (strcmp(*(*(bonds + index) + i), ""))
      {
        bondinds[bondi] = i;
        bondtypes[bondi] = *(*(bonds + index) + i);
        bondi++;
      }
    }
    int maxleaves = (int)pow((double)MAXBONDS, (double)depth); // Maximum possible number of leaf nodes at this depth
    char **outlist = malloc(sizeof(char *) * (maxleaves + 1));
    if (!outlist)
    { // If the outlist pointer wasn't allocated, you've hit the recursion limit
      return NULL;
    }
    *outlist = NULL;
    int leafind = 0;
    for (i = 0; i < bondi; i++)
    {
      if (bondinds[i] != prevind)
      { // Don't analyze the atom we just came from in the parent function call
        char newpre[strlen(prestring) + 8];
        strcpy(newpre, prestring);
        strcat(newpre, bondtypes[i]);
        strcat(newpre, *(atoms + bondinds[i]));
        // Recurse and fetch all leaves of the binding tree for this neighbor
        char **new = buildTree(depth - 1, bondinds[i], atoms, bonds, newpre, index, atomcount);
        char **newit = new;
        while (*newit)
        { // Append the new leaves onto outlist
          *(outlist + leafind) = *newit;
          newit++;
          leafind++;
        }
        free(new);
        *(outlist + leafind) = NULL;
      }
    }
    if (!*outlist)
    { // Another base case, if the current atom's only neighbor is the atom analyzed in the parent function call
      free(outlist);
      outlist = malloc(sizeof(char *) * 2);
      *outlist = malloc(sizeof(char) * (strlen(prestring) + 1));
      strcpy(*outlist, prestring);
      *(outlist + 1) = NULL;
    }
    return outlist;
  }
}

double searchAssigns(int atomcount, int **allcands, int candcounts[], int *assign, char ***tempbond, char ***querybond, double **querycoord, double **tempcoord, int *bestassign)
{
  double **dists = malloc(sizeof(double *) * atomcount); // Distances between query atoms and template atoms
  double querydists[atomcount * atomcount];              // Distances between all query atoms
  int queryconnect[atomcount][MAXBONDS];
  int bondcount[atomcount];
  int connectcount[atomcount];
  int index;
  int i;
  int j;

  // precalculate all query-template atomic distances
  for (i = 0; i < atomcount; i++)
  {
    *(assign + i) = -1;
    connectcount[i] = 0;
    double *distind = malloc(sizeof(double) * candcounts[i]);
    for (j = 0; j < candcounts[i]; j++)
    {
      double dist = 0.0;
      for (index = 0; index < 3; index++)
      {
        dist += pow(*(*(querycoord + i) + index) - *(*(tempcoord + *(*(allcands + i) + j)) + index), 2.0);
      }
      *(distind + j) = dist;
    }
    *(dists + i) = distind;
  }
  memcpy(bestassign, assign, sizeof(int) * atomcount);

  // precalculate all query-query atomic distances for feasibility check
  for (i = 0; i < atomcount; i++)
  {
    for (j = i; j < atomcount; j++)
    {
      double dist = 0.0;
      for (index = 0; index < 3; index++)
      {
        dist += pow(*(*(querycoord + i) + index) - *(*(querycoord + j) + index), 2.0);
      }
      querydists[i * atomcount + j] = pow(dist, 0.5);
      querydists[j * atomcount + i] = querydists[i * atomcount + j];
    }
  }

  // Calculate bond degree for every atom
  for (i = 0; i < atomcount; i++)
  {
    int degree = 0;
    for (j = 0; j < atomcount; j++)
    {
      if (strcmp(*(*(querybond + i) + j), ""))
      {
        queryconnect[i][degree] = j;
        degree++;
      }
    }
    bondcount[i] = degree;
  }

  // bubble sort all possible atoms at each position by query-template distance
  for (index = 0; index < atomcount; index++)
  {
    for (i = 0; i < candcounts[index]; i++)
    {
      for (j = 0; j < candcounts[index] - i - 1; j++)
      {
        if (*(*(dists + index) + j) > *(*(dists + index) + j + 1))
        {
          int tempind = *(*(allcands + index) + j);
          double tempdist = *(*(dists + index) + j);
          *(*(allcands + index) + j) = *(*(allcands + index) + j + 1);
          *(*(dists + index) + j) = *(*(dists + index) + j + 1);
          *(*(allcands + index) + j + 1) = tempind;
          *(*(dists + index) + j + 1) = tempdist;
        }
      }
    }
  }
  /*
  for (i=0;i<atomcount;i++){
    for (j=0;j<candcounts[i];j++){
      printf("%d ",*(*(allcands+i)+j));
    }
    printf("\n");
  }
  */
  int history[atomcount];
  int histinds[atomcount];
  index = 0;
  for (i = 0; i < atomcount; i++)
  {
    histinds[i] = 0;
  }

  double runningTotal = 0.0;
  double bestTotal = DBL_MAX;
  while (1)
  { // While not all mappings have been searched
    if (index == atomcount)
    { // If we've reached the end of a mapping and haven't been pruned
      if (runningTotal < bestTotal)
      {
        memcpy(bestassign, assign, sizeof(int) * atomcount);
        bestTotal = runningTotal;
      }
      // printf("%.2f\n",bestTotal);
      index--;
      runningTotal -= *(*(dists + history[index]) + histinds[index] - 1);
      continue;
    }
    if (histinds[index])
    { // Atom to analyze has been picked, we need to change the mapping

      while (index > 0 && histinds[index] == candcounts[history[index]])
      {
        histinds[index] = 0;
        *(assign + history[index]) = -1;
        for (i = 0; i < bondcount[history[index]]; i++)
        {
          connectcount[queryconnect[history[index]][i]]--;
        }
        index--;
        runningTotal -= *(*(dists + history[index]) + histinds[index] - 1);
      }
      if (index == 0 && histinds[0] == candcounts[history[0]])
      { // This occurs when all mappings have been exhausted
        break;
      }
    }
    else
    { // Pick an atom to analyze
      int nextAtom = 0;
      double bestMetric = DBL_MAX;
      for (i = 0; i < atomcount; i++)
      {
        double nodescore = -0.1 * ((double)connectcount[i]) + 1.0 * ((double)candcounts[i]);
        if (nodescore < bestMetric && *(assign + i) == -1)
        {
          nextAtom = i;
          bestMetric = nodescore;
        }
      }
      history[index] = nextAtom;
      for (i = 0; i < bondcount[history[index]]; i++)
      {
        connectcount[queryconnect[history[index]][i]]++;
      }
      // printf("Next atom: %d\n",nextAtom);
    }
    int foundflag = 0;
    for (i = histinds[index]; i < candcounts[history[index]]; i++)
    {

      if (runningTotal + *(*(dists + history[index]) + i) > bestTotal)
      { // Dead end elimination check
        break;
      }

      if (!inArray(*(*(allcands + history[index]) + i), assign, atomcount) && validateBonds(assign, *(*(allcands + history[index]) + i), history[index], querybond, tempbond, querydists, atomcount))
      { // Feasibility check
        foundflag = 1;
        *(assign + history[index]) = *(*(allcands + history[index]) + i);
        histinds[index] = i + 1;
        runningTotal += *(*(dists + history[index]) + i);
        index++;
        break;
      }
    }
    if (!foundflag)
    { // This occurs if none of the remaining possibilities can be mapped
      if (index == 0)
      {
        break;
      }
      else
      {
        histinds[index] = 0;
        *(assign + history[index]) = -1;
        for (i = 0; i < bondcount[history[index]]; i++)
        {
          connectcount[queryconnect[history[index]][i]]--;
        }
        index--;
        runningTotal -= *(*(dists + history[index]) + histinds[index] - 1);
      }
    }
  }
  for (i = 0; i < atomcount; i++)
  {
    free(*(dists + i));
  }
  free(dists);
  if (*bestassign != -1)
  {
    return pow(bestTotal / ((double)atomcount), 0.5);
  }
  else
  {
    return DBL_MAX;
  }
}

// Checks if the assignment of the current atom is feasible
int validateBonds(int *atomassign, int proposedatom, int assignpos, char ***querybond, char ***tempbond, double querydists[], int atomcount)
{
  int i;
  int j;
  int k;

  int queryinds[MAXBONDS];
  int queri = 0;
  for (j = 0; j < atomcount; j++)
  {
    if (strcmp(*(*(querybond + assignpos) + j), ""))
    {
      queryinds[queri] = j;
      queri++;
    }
  }
  for (j = 0; j < queri; j++)
  {
    int assignatom = *(atomassign + queryinds[j]);
    if (assignatom >= 0 && !strcmp(*(*(tempbond + proposedatom) + assignatom), ""))
    {
      return 0;
    }
  }

  return 1;
}

// Returns the lowest RMSD of all possible mappings for query atoms with template indices given the two molecules' bonding network
double assignAtoms(char **tempatom, char ***tempbond, char **queryatom, char ***querybond, double **querycoord, double **tempcoord, int *querynums, int *tempnums, int atomcount, int simpleflag)
{
  int **allcands = malloc(sizeof(int *) * atomcount); // List of all atoms in the template that could feasibly be each query atom
  int candcounts[atomcount];                          // Number of atoms in the template that could feasibly be each query atom
  int i;

  // Iterate through each query atom and determine which template atoms correspond to the query
  for (i = 0; i < atomcount; i++)
  {
    int j;
    int candidates[atomcount]; // Flags corresponding to if each template atom could correspond to the current query atom
    int viablecands = 0;       // Count of template atoms that could correspond to the current query atom
    for (j = 0; j < atomcount; j++)
    {
      if (!strcmp(*(queryatom + i), *(tempatom + j)))
      {
        candidates[j] = 1;
        viablecands++;
      }
      else
      {
        candidates[j] = 0;
      }
    }
    int treedepth = 1; // Recursion depth
    while (treedepth <= MAXDEPTH)
    { // Recurse deeper until you've searched all atoms or you've hit the recursion limit
      char **qtree = buildTree(treedepth, i, queryatom, querybond, *(queryatom + i), -1, atomcount);
      if (!qtree)
      { // This means you've hit the recursion limit and can't analyze any deeper
        break;
      }
      char **qit = qtree;
      int qlen = 0;
      while (*qit)
      {
        qlen++;
        qit++;
      }
      for (j = 0; j < atomcount; j++)
      {
        if (candidates[j])
        {
          char **ttree = buildTree(treedepth, j, tempatom, tempbond, *(tempatom + j), -1, atomcount);
          char **tit = ttree;
          int tlen = 0;
          while (*tit)
          {
            tlen++;
            tit++;
          }
          /*
          int n;
          printf("Candidate %d:",j);
          for(n=0;n<tlen;n++){
            printf("%s,",*(ttree+n));
          }
          printf("\nQuery %d:",i);
          for(n=0;n<qlen;n++){
            printf("%s,",*(qtree+n));
          }
          printf("\n");
          */
          if (tlen != qlen || !arrayIdentity(qtree, ttree, qlen))
          { // If the template atom tree and query atom tree don't have the same leaves, they're not the same atom
            candidates[j] = 0;
            viablecands--;
          }
          tit = ttree;
          while (*tit)
          {
            free(*tit);
            tit++;
          }
          free(ttree);
        }
      }
      qit = qtree;
      while (*qit)
      {
        free(*qit);
        qit++;
      }
      free(qtree);
      treedepth++;
    }
    if (!viablecands)
    { // If there's no possible atom, something went wrong or the two molecules are not identical
      if (!generalizeBonds(querybond, atomcount))
      {
        if (!simpleflag)
        {
          printf("No atoms mappable for atom %d, generalizing bonds...\n", i);
        }
        generalizeBonds(tempbond, atomcount);
        for (j = 0; j < i; j++)
        {
          free(*(allcands + j));
          candcounts[j] = 0;
        }
        i = -1;
        continue;
      }
      else
      {
        printf("Atom assignment failed for atom %d.\n", i);
        return DBL_MAX;
      }
    }
    else
    { // Otherwise, store all possible template atoms for this query atom
      candcounts[i] = viablecands;
      int *atomcands = malloc(sizeof(int) * viablecands);
      int k = 0;
      for (j = 0; j < atomcount; j++)
      {
        if (candidates[j])
        {
          *(atomcands + k) = j;
          k++;
          if (k == viablecands)
          {
            break;
          }
        }
      }
      *(allcands + i) = atomcands;
    }
  }

  double possiblemaps = 1.0;
  for (i = 0; i < atomcount; i++)
  {
    possiblemaps *= candcounts[i];
  }
  
  // Calculate RMSD of all possible mappings given each query atoms possible template atoms and return the minimum
  int *assign = malloc(sizeof(int) * atomcount);
  int *bestassign = malloc(sizeof(int) * atomcount);
  double bestrmsd = searchAssigns(atomcount, allcands, candcounts, assign, tempbond, querybond, querycoord, tempcoord, bestassign);

  if (bestrmsd == DBL_MAX)
  {
    printf("No valid mapping exists\n");
    return DBL_MAX;
  }

  if (!simpleflag)
  {
    printf("Calculated Docking RMSD: %.3f\n\n", bestrmsd);
    printf("Total # of Possible Mappings: %.0f\n", possiblemaps);
    printf("Optimal mapping (First file -> Second file, * indicates correspondence is not one-to-one):\n");
    for (i = 0; i < atomcount; i++)
    {
      printf("%s%3d -> %s%3d ", *(queryatom + i), *(querynums + i), *(tempatom + *(bestassign + i)), *(tempnums + *(bestassign + i)));
      if (*(querynums + i) == *(tempnums + *(bestassign + i)))
      {
        printf("\n");
      }
      else
      {
        printf("*\n");
      }
    }
  }
  

  free(assign);
  free(bestassign);

  return bestrmsd;
}
