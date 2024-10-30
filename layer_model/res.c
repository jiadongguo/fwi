#include <stdio.h>
int main()
{
	FILE *fd1,*fd2,*fd3;
	fd1=fopen("vini","rb");
	fd2=fopen("vinv","rb");
	fd3=fopen("vres","wb");
	float a[200*300],b[200*300],c[200*300];
	fread(a,sizeof(float)*200*300,1,fd1);
	fread(b,sizeof(float)*200*300,1,fd2);
	for(int i=0;i<200*300;i++)
	{
		c[i]=b[i]-a[i];
	}
	fwrite(c,sizeof(float)*200*300,1,fd3);
	fclose(fd1);
	fclose(fd2);
	fclose(fd3);
	return 0;
}
