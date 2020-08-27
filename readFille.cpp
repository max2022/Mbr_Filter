#include <iostream>
#include <string>
#include <bits/stdc++.h>

using namespace std;

struct table_row {
    int id;
    float x;
    float y;
};

vector<int> seletedFeaturesSizes;
vector<int> seletedFeatures;


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

// read Arpan's outfile as input
void readCombinations(const char *file) {
	string myText;
	int b = 0, c;

	ifstream MyReadFile(file);

	while (getline (MyReadFile, myText)) {
		c = 0;
		while ((b = myText.find(",")) != std::string::npos) {
		    seletedFeatures.push_back(stoi(myText.substr(0, b)));
		    myText.erase(0, b + 1);
		    c++;
		} 
	    seletedFeaturesSizes.push_back(++c);	
		seletedFeatures.push_back(stoi(myText));
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

// 	cout << "---" << endl;

// 	for (int i = 0; i < seletedFeaturesSizes.size(); ++i)
// 	{
// 		cout << seletedFeaturesSizes[i] << endl;
// 	}
// }