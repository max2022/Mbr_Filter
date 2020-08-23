#include <iostream>
#include <string>
#include <bits/stdc++.h>

//#define ROWS 44006
//#define ROWS 165458
//small data
//#define ROWS 650

using namespace std;

struct table_row {
    int id;
    float x;
    float y;
};

int ROWS;

// read file data and retunr data as an object array
struct table_row *createArray(const char *fileName) {
	FILE *fp = fopen(fileName, "r");

	if (!fp) {
        printf("Can not open file\n");
        return NULL;
    }
	
	fscanf(fp, "%d", &ROWS);
	cout << "Total Instances: " << ROWS << endl;
	
	struct table_row* table_rows = (struct table_row*)malloc(sizeof(struct table_row) * ROWS);
	int count = 0;

	
	for (count = 0; count < ROWS; ++count) {
		fscanf(fp, "%d, %f, %f", &table_rows[count].id, &table_rows[count].x, &table_rows[count].y);
		//cout << "id is -> "<< table_rows[count].id << endl;
	}
	//cout << "reading done" << endl;
	fclose(fp);
    return table_rows;
}


