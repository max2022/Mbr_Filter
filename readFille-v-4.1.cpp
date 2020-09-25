#include <iostream>
#include <string>
#include <cmath> 
#include <bits/stdc++.h>

using namespace std;

struct table_row {
    int id;
    float x;
    float y;
};

typedef std::bitset<FMAX> thirteenBits;

vector<thirteenBits> seletedFeatures;


// number of rows in the data file
int ROWS;

// grid number of rows
int GRID_ROWS;

// grid number of columns
int GRID_COLS;

// grid origin
float GRID_MIN_X, GRID_MIN_Y;

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
	float x,y;

	fscanf(fp, "%f, %f, %f, %f", &GRID_MIN_X, &GRID_MIN_Y, &x, &y);

	GRID_COLS = ceil((x - GRID_MIN_X)/(DIST * 2)) + 1;
	GRID_ROWS = ceil((y - GRID_MIN_Y)/(DIST * 2)) + 1;
	cout << GRID_ROWS << " " << GRID_COLS << endl;


	fclose(fp);
    return table_rows;
}

// read Arpan's outfile as input
void readCombinations(const char *file) {
	string myText;
	int b = 0, c;
	thirteenBits n;

	ifstream MyReadFile(file);

	while (getline (MyReadFile, myText)) {
		c = 0;
		n.reset();
		while ((b = myText.find(",")) != std::string::npos) {
			n[FMAX - 1 - stoi(myText.substr(0, b))] = 1;
		    // seletedFeatures.push_back(stoi(myText.substr(0, b)));
		    myText.erase(0, b + 1);
		    c++;
		} 
		n[FMAX - 1 - stoi(myText)] = 1;	    
		seletedFeatures.push_back(n);
	}
	
	MyReadFile.close();
}

// int main() {
// 	vector<int> x;
// 	readCombinations("Arpan-input.txt");
// 	for (int i = 0; i < seletedFeatures.size(); ++i)
// 	{
// 		cout << seletedFeatures[i] << endl;
// 	}
// }
