#include <iostream>
#include <string>

// #define ROWS 165459 
#define ROWS 44006

struct table_row {
    int id;
    float x;
    float y;
};

// read file data and retunr data as an object array
struct table_row *createArray(const char *fileName) {
	FILE *fp = fopen(fileName, "r");

	if (!fp) {
        printf("Can not open file\n");
        return NULL;
    }

	char str[1024];
	struct table_row* table_rows = (struct table_row*)malloc(sizeof(struct table_row) * ROWS);
	int count = 0;

	fscanf(fp, "%s", str);

	for (; count < ROWS; ++count) {
		fscanf(fp, "%d, %f, %f", &table_rows[count].id, &table_rows[count].x, &table_rows[count].y);
	}

	fclose(fp);
    return table_rows;
}

// int main() {
	
// 	struct table_row *dat;
// 	dat = createArray("Seattle2012.csv");

// 	for (int i = 0; i < ROWS; ++i)
// 	{
// 		std::cout << dat[i].id << ": " << dat[i].x << ": " << dat[i].y << std::endl;
// 	}

// 	return 0;
// }