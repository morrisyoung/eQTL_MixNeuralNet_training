

// standard libraries:
#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <string.h>
#include <string>
#include <array>
#include <forward_list>
#include <utility>
#include <vector>
#include <sys/time.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <math.h>


using namespace std;



int main()
{

	//=========== write the file ============
	char filename[100] = "try.txt";
    FILE * file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

    for(int i=0; i<2; i++)
    {
    	//
		float number = 1234123412341234.5678901234;
    	fwrite(&number, 1, sizeof(number), file_out);
		fwrite("\t", 1, sizeof(char), file_out);
    	fwrite(&number, 1, sizeof(number), file_out);
		fwrite("\t", 1, sizeof(char), file_out);
    	fwrite(&number, 1, sizeof(number), file_out);
		fwrite("\t", 1, sizeof(char), file_out);
    	fwrite(&number, 1, sizeof(number), file_out);
		fwrite("\n", 1, sizeof(char), file_out);
    }


	//char buf[1024];
	//sprintf(buf, "%.20f\t", number);
	//fwrite(buf, sizeof(char), strlen(buf), file_out);
	//fwrite(&number, sizeof(int), 1, file_out);
	//fprintf(file_out, "%f", number);

	fclose(file_out);


	printf("%.8f\n", 123412341234.5678901234);



	//=========== read the file =============
    FILE * file_in = fopen(filename, "r");
    if(file_in == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }
    for(int i=0; i<2; i++)
    {
    	//
		float number;
		char c;
    	fread(&number, 1, sizeof(float), file_in);
    	fread(&c, 1, sizeof(char), file_in);
    	cout << number << endl;
       	fread(&number, 1, sizeof(float), file_in);
    	fread(&c, 1, sizeof(char), file_in);
    	cout << number << endl;
       	fread(&number, 1, sizeof(float), file_in);
    	fread(&c, 1, sizeof(char), file_in);
    	cout << number << endl;
       	fread(&number, 1, sizeof(float), file_in);
    	fread(&c, 1, sizeof(char), file_in);
    	cout << number << endl;
    	cout << "**" << endl;
    }
	fclose(file_in);




	//================ other stuff ================
	char array[100] = "nan";
	float number = stof(array);
	cout << number << endl;

	cout << isnan(number) << endl;
	number = 1.2;
	cout << isnan(number) << endl;
	cout << number << endl;
	cout << isnan(12.34) << endl;






	return 0;
}

