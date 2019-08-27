/****************************************
 Parallelization of Principal Component Analysis
 By Kunal Jani
*****************************************/
#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<string.h>
#include<math.h>
#include<time.h>
#define N 3
/********************
struct eigenvector
{
	char *value;
};
struct eigenvalue
{
	char *value;
};
***********************/
int main()
{
	int rows=24;
	int cols=4;
	double a[rows][cols];
	double b[rows][cols];
	double means[cols];
	int i,j;
	double sum;
	//Assigning values to the original matrix.
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			a[i][j]=i+j;
	//		printf("%f ",a[i][j]);
		}
	//	printf("\n");
	}
	//Calculating mean for each column.
	for(i=0;i<cols;i++)
	{
		sum=0;
		for(j=0;j<rows;j++)
		{
			sum=sum+a[j][i];
		}
		means[i]=sum/rows;
	//	printf("%f\n",means[i]);
	}
	//Suntrracting the mean from each matrix,
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			b[i][j]=a[i][j]-means[j];
	//		printf("%f ",b[i][j]);
		}
	//	printf("\n");
	}
	//Obtaining he covariance matrix and writing it to a file.
	double c[cols][rows];
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			c[j][i]=b[i][j];
		}
	}
	/*
	for(i=0;i<cols;i++)
	{
		for(j=0;j<rows;j++)
		{
			printf("%f ",c[i][j]);
		}
		printf("\n");
	}
	*/
	double d[cols][cols];
	int k;
	FILE *fp;
	fp=fopen("matrix.csv","w+");
	for(i=0;i<cols;i++)
	{
		fprintf(fp,"Element,");
	}
	fprintf(fp,"\n");
	for(i=0;i<cols;i++)
	{
		for(j=0;j<cols;j++)
		{
			d[i][j]=0;
			for(k=0;k<rows;k++)
			{
				d[i][j]=d[i][j]+c[i][k]*b[k][j];
			}
			d[i][j]=d[i][j]/(rows-1);
	//		printf("%f,",d[i][j]);
			if(j!=cols-1)
			{
				fprintf(fp,"%f,",d[i][j]);
			}
			else
			{
				fprintf(fp,"%f",d[i][j]);
			}
		}
		fprintf(fp,"\n");
	//	printf("\n");
	}
	fclose(fp);
	//Execution the python script for the finding the eigenvectors and the eigenvalues of the covaraince matrix.
	system("python eigs.py");
	double eigenvectors[cols][cols];
	double eigenvalues[cols];
	char *number;
	size_t len;
	FILE *f1,*f2;
	f1=fopen("eigenvectors.txt","r");
	f2=fopen("eigenvalues.txt","r");
	//Reading the eigenvectors and eigenvalues calculated by the python script.
	for(i=0;i<cols;i++)
	{
		for(j=0;j<cols;j++)
		{
			getline(&number,&len,f1);
	//		printf("%s\n",number);
			eigenvectors[i][j]=(double)atof(number);
	//		printf("%f\n",eigenvectors[i][j]);
		}
	}
	for(i=0;i<cols;i++)
	{
		getline(&number,&len,f2);
	//	printf("%s\n",number);
		eigenvalues[i]=(double)atof(number);
	//	printf("%f\n",eigenvalues[i]);
	}
	fclose(f1);
	fclose(f2);
	//double max1=-30000;
	//Using bubble sort to sort the eigenvectors and the eigenvalues from maximum to minimum.
	for(i=0;i<cols;i++)
	{
		for(j=0;j<cols-1;j++)
		{
			if(eigenvalues[j]<eigenvalues[j+1])
			{
				double tmp=eigenvalues[j];
				eigenvalues[j]=eigenvalues[j+1];
				eigenvalues[j+1]=tmp;
				double tmp2[cols];
				for(k=0;k<cols;k++)
				{
					tmp2[k]=eigenvectors[j][k];
					eigenvectors[j][k]=eigenvectors[j+1][k];
					eigenvectors[j+1][k]=tmp2[k];
				}
			}
		}
	}
/*************
	for(i=0;i<cols;i++)
	{
		printf("%f\n ",eigenvalues[i]);
	}
	

	for(i=0;i<cols;i++)
	{
		for(j=0;j<cols;j++)
		{
			printf("%f ",eigenvectors[i][j]);
		}
		printf("\n");
	}
*************/
	/******
	for(i=0;i<cols;i++)
	{
		for(j=0;j<rows;j++)
		{
			printf("%f ",c[i][j]);
		}
		printf("\n");
	}
	**********/
	double eigenvectors_transpose[cols][cols];
	for(i=0;i<cols;i++)
	{
		for(j=0;j<cols;j++)
		{
			eigenvectors_transpose[j][i]=eigenvectors[i][j];
		}
	}

	double resultant[rows][cols];
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			resultant[i][j]=0;
		}
	}
	
	//Resultant matrix in 2 dimensions obtained by multiplying the first 2 eigenvector rows with
	//the transpose of the matrix obtained by subtracting the mean for each row.

	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			for(k=0;k<cols;k++)
			{
				resultant[i][j]=resultant[i][j]+b[i][k]*eigenvectors[k][j];
			}
		}
	}

/***********
	printf("Resultant matrix first row consists of the x coordinates and the second row consists of ycoordinates to be plotte in python.");
	for(i=0;i<rows;i++)
	{
		for(j=0;j<4;j++)
		{
			printf("%f ",resultant[i][j]);
		}
		printf("\n");
	}
************/
	double resultant2[rows][cols];
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			resultant2[i][j]=0;
		}
	}
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			for(k=0;k<cols;k++)
			{
				resultant2[i][j]=resultant2[i][j]+resultant[i][k]*eigenvectors_transpose[k][j];
			}
		}
	}

//	printf("Resultant matrix first row consists of the x coordinates and the second row consists of ycoordinates to be plotte in python.");
/*********	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			printf("%f ",resultant2[i][j]);
		}
		printf("\n");
	}
***********/
//	printf("\n");

	double new_data[rows][cols];
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			new_data[i][j]=resultant2[i][j]+means[j];
		}
	}

	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			printf("%f ",new_data[i][j]);
		}
		printf("\n");
	}


/**********
	double eigenvectors_transposed[cols][2];
	for(i=0;i<cols;i++)
	{
		for(j=0;j<2;j++)
		{
			eigenvectors_transposed[j][i]=eigenvectors[i][j];
		}
	}
	double original[cols][rows];
	for(i=0;i<cols;i++)
	{
		for(j=0;j<rows;j++)
		{
			original[i][j]=0;
			for(k=0;k<2;k++)
			{
				original[i][j]=original[i][j]+(eigenvectors_transposed[i][k]*resultant[k][j]);
			}
			printf("%f ",original[i][j]);
		}
		printf("\n");
	}
***********/
/**********
	FILE *fpps;
	fpps=fopen("principal_components.csv","w+");
	fprintf(fpps,",x,y");
	for(i=0;i<rows;i++)
	{
		fprintf(fpps,"\n%d",i+1);
		fprintf(fpps,",%f",resultant[0][i]);
		fprintf(fpps,",%f",resultant[1][i]);
	}
	fclose(fpps);
************/

	FILE *fpps;
	fpps=fopen("principal_components2.csv","w+");
	for(i=0;i<cols;i++)
	{
		fprintf(fpps,",%d",i);
	}
//	fprintf(fpps,",w,z,x,y");
	for(i=0;i<rows;i++)
	{
		fprintf(fpps,"\n%d",i+1);
		for(j=0;j<cols;j++)
		{
			fprintf(fpps,",%f",resultant[j][i]);
		}
	}
	fclose(fpps);	
	return 0;

}
