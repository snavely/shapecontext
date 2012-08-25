/* filter.cpp */
/* Learn pairwise distributions between sets of points */

#include <stdio.h>
#include <stdlib.h>

#define NUM_FEATURES 13

int learnDistrs(int argc, char **argv) 
{
    char *imageList_in;

    if (argc != 2) {
	printf("Usage: %s learnDistrs <imageList.txt>\n", name);
	return -1;
    }

    imageList_in = argv[1];

    /* Read the feature names */
    FILE *f = fopen(imageList_in, "r");
    std::vector<char *> imageList;

    while (fgets(buf, 1024, f) != NULL) {
	imageList.push_back(strdup(buf));
    }

    fclose(f);

    /* Read the feature files */
    typedef std::vector<point_t> point_list_t;
    std::vector<point_list_t> pts;
    
    for (int i = 0; i < (int) imageList.size(); i++) {
	f = fopen(imageList[i], "r");

	/* Read boundaries */
	double xMin, yMin, xMax, yMax;
	fscanf(f, "%f %f %f %f", &xMin, &yMin, &xMax, &yMax);

	point_list_t ptList;

	/* Read the feature positions */
	for (int j = 0; j < NUM_FEATURES; j++) {
	    double x, y;
	    point_t pt;
	    fscanf(f, "%f %f", x, y);

	    pt.x = x;
	    pt.y = y;
	    
	    ptList.push_back(pt);
	}

	pts.push_back(ptList);

	fclose(f);
    }

    /* Compute means and covariances */
    double *means = new double[NUM_FEATURES * 2];
    double *covariances = new double[NUM_FEATURES * NUM_FEATURES * 4];




}

void printUsage(char *name) 
{
    printf("Usage: %s learnDistrs <imageList.txt>\n", name);
    printf("       %s filterPoints <pointSet.txt>\n", name);
}


int main(int argc, char **argv) 
{
    if (argc < 2) {
	printUsage(argv[0]);
	return -1;
    }
    
    if (strcmp(argv[1], "learnDistrs") == 0) {
	return learnDistrs(argc, argv);
    } else if (strcmp(argv[2], "filterPoints") == 0) {
	return filterPoints(argc, argv);
    } else {
	printUsage(argv[0]);
	return -1;
    }

    return 0;
}
