#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<time.h>
#include<string.h>
int main(int argc,char **argv)
{
	int rows=24;
	int cols=4;
	double a[rows][cols];
	double b[rows][cols];
	double c[cols][rows];
	double means[rows];
	int i,j,k;
	double sum;
	int ierr,p,id;
	MPI_Status status;
	ierr=MPI_Init(&argc,&argv);
	ierr=MPI_Comm_rank(MPI_COMM_WORLD,&id);
	ierr=MPI_Comm_size(MPI_COMM_WORLD,&p);
	if(id==0)
	{
		for(i=0;i<rows;i++)
		{	
			for(j=0;j<cols;j++)
			{
				a[i][j]=i+j;
		//		printf("%f ",a[i][j]);
			}
	//	printf("\n");
		}

		for(i=0;i<cols;i++)
		{
			sum=0;
			double elements[rows/(p-1)];
			int cnt=0;
			for(j=1;j<p;j++)
			{
				for(k=0;k<rows/(p-1);k++)
				{
					elements[k]=a[cnt++][i];
				}
				ierr=MPI_Send(elements,rows/(p-1),MPI_DOUBLE,j,0,MPI_COMM_WORLD);
			}
		}
		double sm[cols];
		double summs[cols];
		for(i=0;i<cols;i++)
		{
			summs[i]=0;
		}
		for(j=1;j<p;j++)
		{
			ierr=MPI_Recv(&sm,cols,MPI_DOUBLE,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<cols;i++)
			{
				summs[i]=summs[i]+sm[i];
			}
		}
		
	//	printf("Mean:\n");
		for(i=0;i<cols;i++)
		{
			means[i]=summs[i]/rows;
			printf("%f\n",means[i]);
		}
		for(i=1;i<p;i++)
		{
			ierr=MPI_Send(means,cols,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
		}
		for(i=0;i<cols;i++)
		{
			sum=0;
			double elements[rows/(p-1)];
			int cnt=0;
			for(j=1;j<p;j++)
			{
				for(k=0;k<rows/(p-1);k++)
				{
					elements[k]=a[cnt++][i];
				}
				ierr=MPI_Send(elements,rows/(p-1),MPI_DOUBLE,j,0,MPI_COMM_WORLD);
			}
		}
		double bb[rows/(p-1)];
		int cnt=0;
		for(i=0;i<cols;i++)
		{
			cnt=0;
			for(j=1;j<p;j++)
			{
				ierr=MPI_Recv(&bb,rows/(p-1),MPI_DOUBLE,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				for(k=0;k<rows/(p-1);k++)
				{
					b[cnt++][i]=bb[k];
				}
			}
		}
		/**************
		for(i=0;i<rows;i++)
		{
			for(j=0;j<cols;j++)
			{
				printf("%f ",b[i][j]);
			}
			printf("\n");
		}
		*************/
		for(i=0;i<rows;i++)
		{
			for(j=0;j<cols;j++)
			{
				c[j][i]=b[i][j];
			}
		}
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
				for(k=0;k<cols;k++)
				{
					d[i][j]=d[i][j]+c[i][k]*b[k][j];
				}
				d[i][j]=d[i][j]/(rows-1);
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
		}
		fclose(fp);
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
	double resultant[2][rows];
	for(i=0;i<2;i++)
	{
		for(j=0;j<rows;j++)
		{
			resultant[i][j]=0;
		}
	}
	for(i=0;i<2;i++)
	{
		for(j=1;j<p;j++)
		{
			ierr=MPI_Send(eigenvectors[i],cols,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
		}
	}
	cnt=0;
	for(i=1;i<p;i++)
	{
		for(j=0;j<rows/(p-1);j++)
		{
			ierr=MPI_Send(b[cnt++],cols,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
		}
	}
	double result[rows/(p-1)];
	for(i=0;i<2;i++)
	{
		int counter=0;
		for(j=1;j<p;j++)
		{
			ierr=MPI_Recv(&result,rows/(p-1),MPI_DOUBLE,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(k=0;k<rows/(p-1);k++)
			{
				resultant[i][counter++]=result[k];
			}
		}
	}
	for(i=0;i<2;i++)
	{
		for(j=0;j<rows;j++)
		{
			printf("%f ",resultant[i][j]);
		}
		printf("\n");
	}
	}
	else
	{
		double sums[cols];
		for(i=0;i<cols;i++)
		{
			sums[i]=0;
		}
		double els[rows/(p-1)];
		for(i=0;i<cols;i++)
		{
			ierr=MPI_Recv(&els,rows/(p-1),MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(j=0;j<rows/(p-1);j++)
			{
				sums[i]=sums[i]+els[j];
			//	printf("%f ",els[j]);
			}
		//	printf("\n");
		}
		/**************
		printf("Sums:\n");
		
		for(i=0;i<cols;i++)
		{
			printf("%f\n",sums[i]);
		}
		******************/
		ierr=MPI_Send(sums,cols,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		ierr=MPI_Recv(&means,cols,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	/*********	
		printf("Means:\n");
		for(i=0;i<cols;i++)
		{
			printf("%f\n",means[i]);
		}
	************/
		double bb[rows/(p-1)];	
		for(i=0;i<cols;i++)
		{
			ierr=MPI_Recv(&els,rows/(p-1),MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(j=0;j<rows/(p-1);j++)
			{
				bb[j]=els[j]-means[i];
			}
			ierr=MPI_Send(bb,rows/(p-1),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}
		double eigenvectors[2][cols],eigs[cols];
		for(i=0;i<2;i++)
		{
			ierr=MPI_Recv(&eigs,cols,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(j=0;j<cols;j++)
			{
				eigenvectors[i][j]=eigs[j];
			}
		}
		/***********
		for(i=0;i<2;i++)
		{
		//	ierr=MPI_Recv(&eigs,cols,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(j=0;j<cols;j++)
			{
				printf("%f ",eigenvectors[i][j]);
			}
			printf("\n");
		}
		*************/
		double bbb[cols];
		double ccc[rows/(p-1)][cols];
		double result[2][rows/(p-1)];
		for(i=0;i<2;i++)
		{
			for(j=0;j<rows/(p-1);j++)
			{
				result[i][j]=0;
			}
		}
		for(i=0;i<rows/(p-1);i++)
		{
			ierr=MPI_Recv(&bbb,cols,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(j=0;j<cols;j++)
			{
				ccc[i][j]=bbb[j];
		//		printf("%f ",ccc[i][j]);
			}
		//	printf("\n");
		}
		for(i=0;i<2;i++)
		{
			for(j=0;j<rows/(p-1);j++)
			{
				for(k=0;k<cols;k++)
				{
					result[i][j]=result[i][j]+eigenvectors[i][k]*ccc[j][k];
				}
			//	printf("%f ",result[i][j]);
			}
		//	printf("\n");
		}
		for(i=0;i<2;i++)
		{
			ierr=MPI_Send(result[i],rows/(p-1),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}	
	}
	ierr=MPI_Finalize();
}
