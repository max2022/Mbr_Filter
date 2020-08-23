#include <iostream>
#include <string>
#include <bits/stdc++.h>

using namespace std;

struct table_row {
    int id;
    float x;
    float y;
};

// number of rows in the data file
int ROWS;

// read file data and retunr data as an object array
struct table_row *createArray(const char *fileName) {
	FILE *fp = fopen(fileName, "r");

	if (!fp) {
        printf("Can not open file\n");
        return NULL;
    }
	
	fscanf(fp, "%d", &ROWS);
	cout << "Total rows: " << ROWS << endl;

	struct table_row* table_rows = (struct table_row*)malloc(sizeof(struct table_row) * ROWS);
	
	for (int count = 0; count < ROWS; ++count) {
		fscanf(fp, "%d, %f, %f", &table_rows[count].id, &table_rows[count].x, &table_rows[count].y);

	}
	fclose(fp);
    return table_rows;
}

